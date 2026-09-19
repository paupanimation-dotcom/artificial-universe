//! Saving a universe you cannot afford to write down.
//!
//! ARCHITECTURE §19 splits state three ways, and the save format follows it
//! exactly:
//!
//! | Tier | What it is | Cost |
//! |---|---|---|
//! | **Seed** | The entire pristine universe | 8 bytes |
//! | **Dirty chunks** | Everything the simulation actually changed | grows with activity |
//! | **Event log** | The causal record — extinctions, speciations | sparse |
//!
//! A world where a hundred million years have passed but only a few thousand
//! chunks were ever touched saves in kilobytes, because the rest of it is
//! *derivable*. Loading regenerates the pristine parts on demand.
//!
//! The correctness bar is high and it is testable: **a world saved at tick K
//! and resumed must be bit-identical, at tick N, to a world that ran straight
//! through.** If that ever fails, the simulation has hidden state somewhere,
//! and hidden state means the past cannot be trusted. There is a test.

use au_core::event::{Event, EventLog};
use au_core::layer::Layer;
use au_core::time::{SimClock, SimDuration, SimInstant, Tick};

use crate::chunk::{Chunk, ChunkCoord, ChunkStore};
use crate::codec::{CodecError, Reader, Writer, CodecError as CE, FORMAT_VERSION, MAGIC};
use crate::column::ColumnRegistry;

/// Everything needed to resurrect a universe.
pub struct Snapshot {
    pub world_seed: u64,
    pub clock: SimClock,
    pub chunks: Vec<Chunk>,
    pub events: Vec<Event>,
    pub event_cap: usize,
    pub event_total: u64,

    /// Which chunks are under simulation.
    ///
    /// Must be persisted, and it took a moment to see why: it is not derivable
    /// from the chunks themselves. Two worlds with identical contents but
    /// different active sets will evolve differently, because one of them is
    /// simulating a region the other is only storing. Lose this and a reloaded
    /// world quietly stops simulating half of itself.
    pub active: Vec<ChunkCoord>,

    /// The conserved-quantity ledger: what has crossed the boundary of the world.
    /// Four opaque integers — `au-data` neither knows nor needs to know that they
    /// mean energy-in, energy-out, mass-in, mass-out. That interpretation belongs
    /// to the physics layer.
    pub conserved: [i128; 8],
}

fn checksum(bytes: &[u8]) -> u64 {
    let mut h = au_core::hash::Hasher::new();
    h.write_bytes(bytes);
    h.finish()
}

pub fn encode(
    world_seed: u64,
    clock: &SimClock,
    store: &ChunkStore,
    log: &EventLog,
    active: &[ChunkCoord],
    conserved: [i128; 8],
) -> Vec<u8> {
    let mut w = Writer::new();
    w.bytes(MAGIC);
    w.u16(FORMAT_VERSION);
    w.u64(world_seed);

    // Clock: tick, instant, scale — all exact integers, all restorable.
    w.u64(clock.tick().0);
    w.u64(clock.now().secs).u64(clock.now().frac);
    w.u64(clock.scale().secs).u64(clock.scale().frac);

    // Only dirty chunks. This is the whole point.
    let dirty: Vec<&Chunk> = store.dirty_chunks().collect();
    w.u32(dirty.len() as u32);
    for c in dirty {
        c.encode(&mut w);
    }

    // Order: total, capacity, count. Capacity is NOT the count — writing the
    // count here (which I did, first time) restores a log that can hold exactly
    // as many events as it currently has, so it begins throwing away history on
    // the very next emit. The universe resumes and quietly starts forgetting.
    w.u64(log.total());
    w.u64(log.cap() as u64);
    w.u64(log.events().len() as u64);
    for e in log.events() {
        w.u64(e.tick.0);
        w.u8(e.layer as u8);
        w.u16(e.kind);
        w.u64(e.a).u64(e.b).u64(e.c);
    }

    w.u32(active.len() as u32);
    for c in active {
        w.i32(c.x).i32(c.y).i32(c.z);
    }
    for v in conserved {
        w.i128(v);
    }

    let sum = checksum(w.as_slice());
    w.u64(sum);
    w.finish()
}

pub fn decode(bytes: &[u8], reg: &ColumnRegistry) -> Result<Snapshot, CodecError> {
    if bytes.len() < 8 {
        return Err(CE::UnexpectedEof { needed: 8, had: bytes.len() });
    }
    // Verify integrity before trusting a single field.
    let (body, tail) = bytes.split_at(bytes.len() - 8);
    let found = u64::from_le_bytes(tail.try_into().unwrap());
    let expected = checksum(body);
    if found != expected {
        return Err(CE::ChecksumMismatch { expected, found });
    }

    let mut r = Reader::new(body);
    if r.bytes(4)? != MAGIC {
        return Err(CE::BadMagic);
    }
    let version = r.u16()?;
    if version != FORMAT_VERSION {
        // Where migrations will live. Deliberately a hard error today: a
        // best-effort load of an unknown format produces a world that looks
        // right and is wrong, which is worse than not loading at all.
        return Err(CE::UnsupportedVersion(version));
    }

    let world_seed = r.u64()?;
    let tick = Tick(r.u64()?);
    let now = SimInstant { secs: r.u64()?, frac: r.u64()? };
    let scale = SimDuration { secs: r.u64()?, frac: r.u64()? };
    let clock = SimClock::restore(tick, now, scale);

    let n_chunks = r.u32()? as usize;
    let mut chunks = Vec::with_capacity(n_chunks);
    for _ in 0..n_chunks {
        chunks.push(Chunk::decode(&mut r, reg)?);
    }

    let event_total = r.u64()?;
    let event_cap = r.u64()? as usize;
    let n_events = r.u64()? as usize;
    let mut events = Vec::with_capacity(n_events);
    for _ in 0..n_events {
        let tick = Tick(r.u64()?);
        let layer_raw = r.u8()?;
        let layer = *Layer::ALL
            .get(layer_raw as usize)
            .ok_or(CE::Invalid("bad layer id"))?;
        let kind = r.u16()?;
        events.push(Event {
            tick,
            layer,
            kind,
            a: r.u64()?,
            b: r.u64()?,
            c: r.u64()?,
        });
    }

    let n_active = r.u32()? as usize;
    let mut active = Vec::with_capacity(n_active);
    for _ in 0..n_active {
        active.push(ChunkCoord::new(r.i32()?, r.i32()?, r.i32()?));
    }
    let mut conserved = [0i128; 8];
    for c in conserved.iter_mut() {
        *c = r.i128()?;
    }

    Ok(Snapshot { world_seed, clock, chunks, events, event_cap, event_total, active, conserved })
}

impl Snapshot {
    pub fn into_parts(self) -> (u64, SimClock, ChunkStore, EventLog) {
        let mut store = ChunkStore::new(self.world_seed);
        for c in self.chunks {
            store.insert(c);
        }
        let log = EventLog::restore(self.events, self.event_cap.max(1), self.event_total);
        (self.world_seed, self.clock, store, log)
    }
}
