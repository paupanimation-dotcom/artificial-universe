//! Chunks, level of detail, and the memory strategy that makes a planet fit in RAM.
//!
//! ARCHITECTURE §19 and §20 are the load-bearing constraints:
//!
//! > "Avoid storing everything permanently. [...] Procedural Reconstruction:
//! >  generate details when needed. This prevents impossible memory
//! >  requirements."
//!
//! The trick, and it is the whole trick:
//!
//! **A chunk that has never been touched is not data. It is a function of the
//! seed.**
//!
//! So we store only what the simulation has actually *changed*. A pristine
//! chunk can be dropped from memory and regenerated bit-identically later,
//! because generation is `f(world_seed, coord)` and nothing else. A chunk the
//! simulation has written to is dirty, and dirty chunks are the only thing a
//! save file contains.
//!
//! This is what turns "the universe must be enormous" from a wish into an
//! invariant. It only works because [`au_core::Rng`] is derived rather than
//! drawn — regeneration must not depend on *when* it happens.

use au_core::hash::{Hasher, WorldHash};
use au_core::rng::{Domain, Rng};
use std::collections::BTreeMap;

use crate::codec::{CodecError, Reader, Writer};
use crate::column::{ColumnRegistry, ColumnSet};

/// Where a chunk is. Integer lattice; the world is unbounded in principle.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Default)]
pub struct ChunkCoord {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

impl ChunkCoord {
    pub fn new(x: i32, y: i32, z: i32) -> Self {
        ChunkCoord { x, y, z }
    }

    /// Pack into a u64 for event payloads. Lossy on the top bits of each axis;
    /// events are diagnostics, not authoritative storage.
    pub fn packed(self) -> u64 {
        ((self.x as u32 as u64) << 42) | ((self.y as u32 as u64) << 21) | (self.z as u32 as u64 & 0x1F_FFFF)
    }

    /// The RNG key for this chunk. Two axes of key give us (coord, tick) later
    /// without reshaping the API.
    pub fn rng_key(self) -> u64 {
        (self.x as u32 as u64) << 32 | (self.y as u32 as u64) << 16 | (self.z as u32 as u64 & 0xFFFF)
    }
}

/// How finely a chunk is simulated (ARCHITECTURE §20).
///
/// This is not a rendering setting. It selects *which physics is real* in a
/// region. The player can zoom from one organism to a whole planet, and the
/// engine pays for detail only where someone is looking — or where something
/// interesting is happening.
///
/// The hard problem, which we are not solving today but must not design
/// ourselves out of: transitions between levels have to conserve quantities.
/// If demoting a chunk to `Statistical` loses biomass, ecosystems will drift
/// every time the camera moves. Every LOD transition must be a *summary*, never
/// a truncation.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
#[repr(u8)]
pub enum Lod {
    /// Every individual, full physics. Expensive. Rare.
    Full = 0,
    /// Cohorts, not individuals. Simplified interactions.
    Population = 1,
    /// Aggregate mathematics only. Population statistics, no bodies.
    Statistical = 2,
}

impl Lod {
    pub fn from_u8(v: u8) -> Option<Lod> {
        match v {
            0 => Some(Lod::Full),
            1 => Some(Lod::Population),
            2 => Some(Lod::Statistical),
            _ => None,
        }
    }
}

/// One cell of the world.
#[derive(Clone)]
pub struct Chunk {
    pub coord: ChunkCoord,
    pub lod: Lod,
    /// Columns of simulation state. Empty in Phase 1 — there is nothing to
    /// simulate yet, and inventing a payload now would be inventing content.
    pub columns: ColumnSet,
    /// Has the simulation written to this chunk since it was generated?
    ///
    /// Clean ⇒ regenerable from the seed ⇒ evictable, and never saved.
    dirty: bool,
}

impl Chunk {
    pub fn new(coord: ChunkCoord, lod: Lod) -> Self {
        Chunk { coord, lod, columns: ColumnSet::new(), dirty: false }
    }

    #[inline]
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    /// Take a mutable view, marking the chunk as no longer regenerable.
    ///
    /// The only door to mutation. If you can change a chunk without going
    /// through here, eviction will silently destroy the change — and the world
    /// will "forget" things for no visible reason.
    #[inline]
    pub fn columns_mut(&mut self) -> &mut ColumnSet {
        self.dirty = true;
        &mut self.columns
    }

    /// Read-only access. Does not dirty.
    #[inline]
    pub fn columns(&self) -> &ColumnSet {
        &self.columns
    }

    /// LOD is simulation state, so changing it dirties the chunk.
    pub fn set_lod(&mut self, lod: Lod) {
        if self.lod != lod {
            self.lod = lod;
            self.dirty = true;
        }
    }

    pub fn encode(&self, w: &mut Writer) {
        w.i32(self.coord.x).i32(self.coord.y).i32(self.coord.z);
        w.u8(self.lod as u8);
        self.columns.encode(w);
    }

    pub fn decode(r: &mut Reader, reg: &ColumnRegistry) -> Result<Chunk, CodecError> {
        let coord = ChunkCoord::new(r.i32()?, r.i32()?, r.i32()?);
        let lod = Lod::from_u8(r.u8()?).ok_or(CodecError::Invalid("bad LOD"))?;
        let columns = ColumnSet::decode(r, reg)?;
        // A chunk that was saved was dirty by definition — clean ones are not
        // written. Restoring it clean would make it evictable and it would be
        // regenerated back to its pristine state, losing history.
        Ok(Chunk { coord, lod, columns, dirty: true })
    }
}

impl WorldHash for Chunk {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_i32(self.coord.x);
        h.write_i32(self.coord.y);
        h.write_i32(self.coord.z);
        h.write_u8(self.lod as u8);
        self.columns.hash_into(h);
    }
}

/// Builds a chunk from nothing but the seed and the coordinate.
///
/// **Must be a pure function.** No clock, no neighbours, no global state, no
/// "and then ask the ecology system". If generation depends on anything but its
/// inputs, eviction stops being safe and the entire memory strategy collapses.
///
/// Phase 4 will implement the real one (terrain, atmosphere). It will still be
/// pure.
pub trait ChunkGenerator: Send + Sync {
    fn generate(&self, coord: ChunkCoord, rng: &mut Rng) -> Chunk;
}

/// The Phase 1 generator: makes empty chunks.
///
/// This is not a placeholder to be filled with "starter content". The world is
/// genuinely empty because matter does not exist yet — that is Phase 3. The
/// container is real and tested; what it contains is not our decision to make.
pub struct EmptyGenerator;

impl ChunkGenerator for EmptyGenerator {
    fn generate(&self, coord: ChunkCoord, _rng: &mut Rng) -> Chunk {
        Chunk::new(coord, Lod::Statistical)
    }
}

/// The resident set of chunks.
///
/// `BTreeMap` for deterministic iteration — see the note in `column.rs`; the
/// same reasoning applies and the same bug is available if we forget it.
pub struct ChunkStore {
    chunks: BTreeMap<ChunkCoord, Chunk>,
    world_seed: u64,
    /// Diagnostics: how hard is the world working to keep itself in memory?
    pub stat_generated: u64,
    pub stat_evicted: u64,
}

impl ChunkStore {
    pub fn new(world_seed: u64) -> Self {
        ChunkStore {
            chunks: BTreeMap::new(),
            world_seed,
            stat_generated: 0,
            stat_evicted: 0,
        }
    }

    pub fn get(&self, coord: ChunkCoord) -> Option<&Chunk> {
        self.chunks.get(&coord)
    }

    /// Fetch, generating on demand. The normal way to touch the world.
    pub fn get_or_generate(
        &mut self,
        coord: ChunkCoord,
        gen: &dyn ChunkGenerator,
    ) -> &mut Chunk {
        if !self.chunks.contains_key(&coord) {
            let mut rng = Rng::derive(self.world_seed, Domain::CHUNK_GEN, coord.rng_key(), 0);
            let chunk = gen.generate(coord, &mut rng);
            debug_assert!(
                !chunk.is_dirty(),
                "a generator must not produce a dirty chunk; \
                 generation is supposed to be reproducible from the seed"
            );
            self.chunks.insert(coord, chunk);
            self.stat_generated += 1;
        }
        self.chunks.get_mut(&coord).unwrap()
    }

    pub fn get_mut(&mut self, coord: ChunkCoord) -> Option<&mut Chunk> {
        self.chunks.get_mut(&coord)
    }

    /// Drop every clean chunk. Free, lossless, and repeatable.
    ///
    /// This is the memory strategy in one line. If this ever loses information,
    /// something upstream mutated a chunk without going through `columns_mut()`.
    pub fn evict_clean(&mut self) -> usize {
        let before = self.chunks.len();
        self.chunks.retain(|_, c| c.is_dirty());
        let n = before - self.chunks.len();
        self.stat_evicted += n as u64;
        n
    }

    pub fn resident(&self) -> usize {
        self.chunks.len()
    }

    pub fn dirty_count(&self) -> usize {
        self.chunks.values().filter(|c| c.is_dirty()).count()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&ChunkCoord, &Chunk)> {
        self.chunks.iter()
    }

    /// Only dirty chunks are persisted. Everything else is a function of the seed.
    pub fn dirty_chunks(&self) -> impl Iterator<Item = &Chunk> {
        self.chunks.values().filter(|c| c.is_dirty())
    }

    pub fn insert(&mut self, chunk: Chunk) {
        self.chunks.insert(chunk.coord, chunk);
    }
}

impl WorldHash for ChunkStore {
    /// Note what is hashed and what is not: **only dirty chunks**.
    ///
    /// This is deliberate and it is the strongest test in the project. Two
    /// worlds are identical if their *changes* are identical — whether or not
    /// they happen to have the same pristine chunks paged in. A run that evicts
    /// aggressively and one that never evicts must hash the same. If they don't,
    /// procedural reconstruction is broken, and we find out in milliseconds
    /// instead of in year 4,000,000.
    fn hash_into(&self, h: &mut Hasher) {
        let dirty: Vec<&Chunk> = self.dirty_chunks().collect();
        h.write_u64(dirty.len() as u64);
        for c in dirty {
            c.hash_into(h);
        }
    }
}
