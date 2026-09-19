//! Protocells — the first individuals.
//!
//! # What changes here
//!
//! Until now every population in this engine has been a *field*: so many
//! molecules of species 7 in grid cell 412. Fields cannot have descendants. They
//! have no edges, so there is nothing for a boundary to enclose; no identity, so
//! nothing to trace; and no parent, so nothing to inherit from. Phase 5h showed
//! that a barrier is worth having in a flow, and Phase 5i step 1 showed when a
//! bag has enough surface to become two bags. Neither could point at *a* bag.
//!
//! A protocell is a bag with a name. It holds its own molecules, exchanges with
//! the medium through its own membrane, and when geometry allows it, becomes two
//! bags that remember which one it was.
//!
//! # Identity, and the promise Phase 1 made
//!
//! `EntityId` is index + generation, and the comment written beside it in Phase 1
//! said why: organisms die, their slots get reused, and if identity were a bare
//! index then a family tree would eventually point at a stranger. That was
//! written before there was anything to allocate. This is the module that
//! finally allocates.
//!
//! A division **frees the parent's id and allocates two fresh ones**, each
//! recording the parent it came from. The parent does not survive as one of its
//! daughters. That is a choice, and the reason is that a lineage in which the
//! parent continues is ambiguous about which daughter *is* the parent — real
//! bacteria have this problem and biologists argue about it. Here the graph is a
//! clean binary tree: every protocell has exactly one parent and at most two
//! children, and no id ever means two different things.
//!
//! # Where the matter is
//!
//! **Species populations now live in two places.** The grid columns hold the
//! bulk medium; protocells hold what they have taken up. Any accounting of what
//! the world contains must sum both, and every conservation test in this phase
//! does. This is the first time in the project that "how much is there" has
//! required looking in more than one place, and it is worth flagging loudly
//! because a check that forgets the protocells will report matter vanishing
//! every time a vesicle buds.

use au_core::hash::{Hasher, WorldHash};
use au_core::ids::EntityId;
use au_data::chunk::ChunkCoord;

/// One bag, and everything that is true about it.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Protocell {
    pub id: EntityId,
    /// The bag this one split from. `EntityId::NONE` for a vesicle that
    /// nucleated out of the medium rather than being born.
    pub parent: EntityId,
    /// The tick it came into existence. With `parent`, this is a family tree.
    pub born_tick: u64,
    pub chunk: ChunkCoord,
    /// Index of the grid cell it floats in.
    pub cell: u32,
    /// What it holds: `(species, count)`, sorted by species, sparse. Sorted
    /// because iteration order must not depend on how it was built.
    pub contents: Vec<(u32, i128)>,
}

impl Protocell {
    pub fn count(&self, species: u32) -> i128 {
        match self.contents.binary_search_by_key(&species, |(s, _)| *s) {
            Ok(i) => self.contents[i].1,
            Err(_) => 0,
        }
    }

    pub fn add(&mut self, species: u32, n: i128) {
        if n == 0 {
            return;
        }
        match self.contents.binary_search_by_key(&species, |(s, _)| *s) {
            Ok(i) => {
                self.contents[i].1 += n;
                if self.contents[i].1 == 0 {
                    self.contents.remove(i);
                }
            }
            Err(i) => self.contents.insert(i, (species, n)),
        }
    }

    /// Total molecules held, of every kind.
    pub fn total(&self) -> i128 {
        self.contents.iter().map(|(_, n)| *n).sum()
    }
}

impl WorldHash for Protocell {
    fn hash_into(&self, h: &mut Hasher) {
        self.id.hash_into(h);
        self.parent.hash_into(h);
        h.write_u64(self.born_tick);
        h.write_i32(self.chunk.x);
        h.write_i32(self.chunk.y);
        h.write_i32(self.chunk.z);
        h.write_u32(self.cell);
        h.write_u32(self.contents.len() as u32);
        for (s, n) in &self.contents {
            h.write_u32(*s);
            h.write_bytes(&n.to_le_bytes());
        }
    }
}

/// Every protocell in the world, in a deterministic order.
///
/// A flat vector sorted by id rather than a map keyed by chunk. Protocells move
/// between chunks the moment anything advects them, and a structure that has to
/// be rebuilt on every move is a structure that will eventually be rebuilt in
/// the wrong order. Sorted by id, iteration is the same everywhere and forever.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct ProtocellStore {
    cells: Vec<Protocell>,
    /// Lifetime counters. History, not state — but hashed, because two worlds
    /// that have seen different numbers of divisions are different worlds even
    /// if they happen to look alike right now.
    pub stat_nucleated: u64,
    pub stat_divided: u64,
    /// Bags that burst — osmotic failure, `v` past the lysis bound.
    pub stat_lysed: u64,
    /// Bags that left through an open port.
    ///
    /// **Kept apart from lysis deliberately.** They are opposite events: one is
    /// a bag failing, the other is a bag being removed without having failed,
    /// and selection is precisely about the difference. Summing them (as this
    /// did when the drain was added) makes a churning population and a turning
    /// over one indistinguishable in the only counters that report either.
    pub stat_washed: u64,
}

impl ProtocellStore {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.cells.len()
    }

    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Protocell> {
        self.cells.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Protocell> {
        self.cells.iter_mut()
    }

    pub fn get(&self, id: EntityId) -> Option<&Protocell> {
        self.cells.binary_search_by_key(&(id.index, id.generation), key).ok().map(|i| &self.cells[i])
    }

    /// Insert, keeping the vector sorted by id.
    pub fn insert(&mut self, p: Protocell) {
        let k = (p.id.index, p.id.generation);
        match self.cells.binary_search_by_key(&k, key) {
            Ok(i) => self.cells[i] = p,
            Err(i) => self.cells.insert(i, p),
        }
    }

    pub fn remove(&mut self, id: EntityId) -> Option<Protocell> {
        let k = (id.index, id.generation);
        match self.cells.binary_search_by_key(&k, key) {
            Ok(i) => Some(self.cells.remove(i)),
            Err(_) => None,
        }
    }

    /// Total molecules of each species held inside protocells, added into `out`.
    ///
    /// The other half of every conservation check in this phase; the grid
    /// columns are the first half.
    pub fn totals_into(&self, out: &mut [i128]) {
        for p in &self.cells {
            for (s, n) in &p.contents {
                if let Some(slot) = out.get_mut(*s as usize) {
                    *slot += *n;
                }
            }
        }
    }

    /// Ancestors of `id`, nearest first, walking `parent` links until they run
    /// out. The family tree the vision promised the player, at the only scale
    /// that currently exists.
    pub fn lineage(&self, id: EntityId) -> Vec<EntityId> {
        let mut out = Vec::new();
        let mut cur = self.get(id).map(|p| p.parent).unwrap_or(EntityId::NONE);
        // A cycle is impossible by construction (a parent always predates its
        // child) but the bound is cheap and a hang is not.
        while !cur.is_none() && out.len() < 4096 {
            out.push(cur);
            match self.get(cur) {
                Some(p) => cur = p.parent,
                None => break, // the ancestor is dead; the trail ends here
            }
        }
        out
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut w = Vec::new();
        w.extend_from_slice(&(self.cells.len() as u64).to_le_bytes());
        for p in &self.cells {
            w.extend_from_slice(&p.id.index.to_le_bytes());
            w.extend_from_slice(&p.id.generation.to_le_bytes());
            w.extend_from_slice(&p.parent.index.to_le_bytes());
            w.extend_from_slice(&p.parent.generation.to_le_bytes());
            w.extend_from_slice(&p.born_tick.to_le_bytes());
            w.extend_from_slice(&p.chunk.x.to_le_bytes());
            w.extend_from_slice(&p.chunk.y.to_le_bytes());
            w.extend_from_slice(&p.chunk.z.to_le_bytes());
            w.extend_from_slice(&p.cell.to_le_bytes());
            w.extend_from_slice(&(p.contents.len() as u32).to_le_bytes());
            for (s, n) in &p.contents {
                w.extend_from_slice(&s.to_le_bytes());
                w.extend_from_slice(&n.to_le_bytes());
            }
        }
        for v in [self.stat_nucleated, self.stat_divided, self.stat_lysed, self.stat_washed] {
            w.extend_from_slice(&v.to_le_bytes());
        }
        w
    }

    pub fn from_bytes(b: &[u8]) -> Option<ProtocellStore> {
        let mut o = 0usize;
        macro_rules! take {
            ($n:expr) => {{
                if o + $n > b.len() {
                    return None;
                }
                let s = &b[o..o + $n];
                o += $n;
                s
            }};
        }
        let n = u64::from_le_bytes(take!(8).try_into().ok()?) as usize;
        let mut cells = Vec::with_capacity(n.min(1 << 20));
        for _ in 0..n {
            let id = EntityId {
                index: u32::from_le_bytes(take!(4).try_into().ok()?),
                generation: u32::from_le_bytes(take!(4).try_into().ok()?),
            };
            let parent = EntityId {
                index: u32::from_le_bytes(take!(4).try_into().ok()?),
                generation: u32::from_le_bytes(take!(4).try_into().ok()?),
            };
            let born_tick = u64::from_le_bytes(take!(8).try_into().ok()?);
            let x = i32::from_le_bytes(take!(4).try_into().ok()?);
            let y = i32::from_le_bytes(take!(4).try_into().ok()?);
            let z = i32::from_le_bytes(take!(4).try_into().ok()?);
            let cell = u32::from_le_bytes(take!(4).try_into().ok()?);
            let k = u32::from_le_bytes(take!(4).try_into().ok()?) as usize;
            let mut contents = Vec::with_capacity(k.min(1 << 16));
            for _ in 0..k {
                let s = u32::from_le_bytes(take!(4).try_into().ok()?);
                let v = i128::from_le_bytes(take!(16).try_into().ok()?);
                contents.push((s, v));
            }
            cells.push(Protocell {
                id,
                parent,
                born_tick,
                chunk: ChunkCoord::new(x, y, z),
                cell,
                contents,
            });
        }
        let stat_nucleated = u64::from_le_bytes(take!(8).try_into().ok()?);
        let stat_divided = u64::from_le_bytes(take!(8).try_into().ok()?);
        let stat_lysed = u64::from_le_bytes(take!(8).try_into().ok()?);
        // Optional: a save written before washout existed has no such field and
        // reads as zero, which is the correct history for a world where no bag
        // could leave.
        let stat_washed = if o + 8 <= b.len() {
            u64::from_le_bytes(take!(8).try_into().ok()?)
        } else {
            0
        };
        let _ = o;
        Some(ProtocellStore { cells, stat_nucleated, stat_divided, stat_lysed, stat_washed })
    }
}

fn key(p: &Protocell) -> (u32, u32) {
    (p.id.index, p.id.generation)
}

impl WorldHash for ProtocellStore {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u64(self.cells.len() as u64);
        for p in &self.cells {
            p.hash_into(h);
        }
        h.write_u64(self.stat_nucleated);
        h.write_u64(self.stat_divided);
        h.write_u64(self.stat_lysed);
        h.write_u64(self.stat_washed);
    }
}
