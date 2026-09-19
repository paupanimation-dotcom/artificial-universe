//! The memory strategy, tested.
//!
//! ARCHITECTURE §19 makes a bold claim: a universe too large to store can be
//! kept in memory anyway, because most of it is a *function of the seed* rather
//! than data. These tests are what stop that from being wishful thinking.
//!
//! The claim has exactly one dangerous failure mode: a chunk gets modified
//! without being marked dirty, gets evicted, gets regenerated in its pristine
//! state — and the world silently forgets something. No crash, no error, just a
//! universe that quietly loses its own history. These tests exist to make that
//! impossible.

use au_core::event::{Event, EventLog};
use au_core::hash::WorldHash;
use au_core::layer::Layer;
use au_core::rng::{Domain, Rng};
use au_core::time::{SimClock, SimDuration, Tick};
use au_data::chunk::*;
use au_data::codec::*;
use au_data::column::*;
use au_data::snapshot;

/// Stand-in for the physics columns Phase 2 will bring. Deliberately generic:
/// we are testing the *container*, and inventing a payload with meaning would
/// be inventing content.
const TEST_FIELD: ColumnId = ColumnId(0xF001);

fn registry() -> ColumnRegistry {
    let mut r = ColumnRegistry::new();
    r.register::<f64>(TEST_FIELD);
    r
}

/// A generator that puts *something* deterministic in each chunk, so that
/// "regenerated identically" is a claim with teeth.
struct NoiseGen;
impl ChunkGenerator for NoiseGen {
    fn generate(&self, coord: ChunkCoord, rng: &mut Rng) -> Chunk {
        let mut c = Chunk::new(coord, Lod::Statistical);
        let mut col = VecColumn::<f64>::new(TEST_FIELD);
        for _ in 0..16 {
            col.push(rng.next_f64());
        }
        // Insert directly, bypassing `columns_mut()`: generated content is NOT
        // a modification. It is derivable, so the chunk stays clean and
        // evictable. This is the one place that distinction matters.
        c.columns.insert(Box::new(col));
        c
    }
}

// ─── Reconstruction ──────────────────────────────────────────────────────────

/// The load-bearing test of the entire memory architecture.
///
/// Generate a chunk. Throw it away. Generate it again. It must be *identical* —
/// not similar, identical — or the universe is not reproducible and "rewind
/// history" is a lie.
#[test]
fn evicted_chunks_regenerate_bit_identically() {
    let mut store = ChunkStore::new(12345);
    let gen = NoiseGen;
    let coords: Vec<_> = (0..50).map(|i| ChunkCoord::new(i % 7, 0, i / 7)).collect();

    for &c in &coords {
        store.get_or_generate(c, &gen);
    }
    let before = store.world_hash();
    let sample: Vec<f64> = store
        .get(coords[9])
        .unwrap()
        .columns()
        .get::<f64>(TEST_FIELD)
        .unwrap()
        .as_slice()
        .to_vec();

    assert_eq!(store.evict_clean(), 50, "nothing was modified, so everything is disposable");
    assert_eq!(store.resident(), 0);

    // Regenerate in a *different order*, to catch any hidden dependence on
    // sequence. There must be none: generation is f(seed, coord).
    for &c in coords.iter().rev() {
        store.get_or_generate(c, &gen);
    }

    assert_eq!(store.world_hash(), before, "the world came back different");
    let after: Vec<f64> = store
        .get(coords[9])
        .unwrap()
        .columns()
        .get::<f64>(TEST_FIELD)
        .unwrap()
        .as_slice()
        .to_vec();
    assert_eq!(sample, after, "a regenerated chunk must be bit-identical");
}

/// The dangerous failure mode, made impossible: anything the simulation touched
/// must survive eviction.
#[test]
fn modified_chunks_are_never_evicted() {
    let mut store = ChunkStore::new(1);
    let gen = NoiseGen;
    for i in 0..10 {
        store.get_or_generate(ChunkCoord::new(i, 0, 0), &gen);
    }

    // Simulate something happening in chunk 3. `columns_mut()` is the only door
    // to mutation, and it dirties.
    {
        let c = store.get_or_generate(ChunkCoord::new(3, 0, 0), &gen);
        assert!(!c.is_dirty(), "pristine until touched");
        c.columns_mut().get_mut::<f64>(TEST_FIELD).unwrap().as_mut_slice()[0] = 999.0;
        assert!(c.is_dirty(), "touching must dirty");
    }

    assert_eq!(store.evict_clean(), 9);
    assert_eq!(store.resident(), 1, "history survives; derivable state does not");
    assert_eq!(
        store.get(ChunkCoord::new(3, 0, 0)).unwrap().columns().get::<f64>(TEST_FIELD).unwrap().as_slice()[0],
        999.0
    );
}

/// The world's identity must not depend on how much RAM the host had.
#[test]
fn world_hash_ignores_memory_pressure() {
    let gen = NoiseGen;

    let mut generous = ChunkStore::new(77);
    let mut cramped = ChunkStore::new(77);

    for i in 0..30 {
        let c = ChunkCoord::new(i, 0, 0);
        generous.get_or_generate(c, &gen);

        cramped.get_or_generate(c, &gen);
        cramped.evict_clean(); // pretend we are desperately short of memory
    }

    // Both worlds changed exactly one chunk, in the same way.
    for store in [&mut generous, &mut cramped] {
        let c = store.get_or_generate(ChunkCoord::new(5, 0, 0), &gen);
        c.columns_mut().get_mut::<f64>(TEST_FIELD).unwrap().as_mut_slice()[0] = 1.0;
    }

    assert_eq!(generous.resident(), 30);
    assert_eq!(cramped.resident(), 1);
    assert_eq!(
        generous.world_hash(),
        cramped.world_hash(),
        "a machine with less RAM must not compute a different universe"
    );
}

// ─── Columns ─────────────────────────────────────────────────────────────────

#[test]
fn columns_round_trip_bit_exactly() {
    // Awkward values on purpose: if the codec ever goes via a decimal string,
    // these are what will catch it.
    let values = vec![0.1, -0.0, 1.0 / 3.0, f64::MIN_POSITIVE, 1e308, -1e-308];
    let mut col = VecColumn::<f64>::new(TEST_FIELD);
    for &v in &values {
        col.push(v);
    }
    let mut set = ColumnSet::new();
    set.insert(Box::new(col));

    let mut w = Writer::new();
    set.encode(&mut w);
    let bytes = w.finish();

    let mut r = Reader::new(&bytes);
    let back = ColumnSet::decode(&mut r, &registry()).unwrap();
    let got = back.get::<f64>(TEST_FIELD).unwrap().as_slice();

    for (a, b) in values.iter().zip(got) {
        assert_eq!(a.to_bits(), b.to_bits(), "float survived only approximately");
    }
    assert_eq!(set.world_hash(), back.world_hash());
}

#[test]
fn an_unregistered_column_is_an_error_not_a_shrug() {
    let mut col = VecColumn::<f64>::new(ColumnId(0xDEAD));
    col.push(1.0);
    let mut set = ColumnSet::new();
    set.insert(Box::new(col));
    let mut w = Writer::new();
    set.encode(&mut w);
    let bytes = w.finish();

    // Silently dropping the column would load a world that *looks* fine and is
    // wrong. That is the worst outcome available, so it must be loud.
    match ColumnSet::decode(&mut Reader::new(&bytes), &registry()) {
        Err(CodecError::UnknownColumn(0xDEAD)) => {}
        Err(e) => panic!("wrong error: {}", e),
        Ok(_) => panic!("an unknown column loaded silently — this is the bad one"),
    }
}

// ─── Snapshots ───────────────────────────────────────────────────────────────

#[test]
fn snapshots_store_changes_not_worlds() {
    let gen = NoiseGen;
    let mut store = ChunkStore::new(5);
    for i in 0..1000 {
        store.get_or_generate(ChunkCoord::new(i % 32, 0, i / 32), &gen);
    }
    // One chunk out of a thousand was actually touched.
    store
        .get_or_generate(ChunkCoord::new(1, 0, 1), &gen)
        .columns_mut()
        .get_mut::<f64>(TEST_FIELD)
        .unwrap()
        .as_mut_slice()[0] = 42.0;

    let clock = SimClock::new(SimDuration::from_secs(1));
    let bytes = snapshot::encode(5, &clock, &store, &EventLog::new(100), &[], [0; 8]);

    // 1000 chunks × 16 f64 = 128 KB of world state. The save should be a rounding
    // error next to that, because 999 of those chunks are derivable.
    assert!(
        bytes.len() < 1000,
        "save was {} bytes; it should hold one dirty chunk, not a thousand pristine ones",
        bytes.len()
    );

    let snap = snapshot::decode(&bytes, &registry()).unwrap();
    assert_eq!(snap.chunks.len(), 1);
    assert_eq!(
        snap.chunks[0].columns().get::<f64>(TEST_FIELD).unwrap().as_slice()[0],
        42.0
    );
}

#[test]
fn snapshots_preserve_clock_and_history() {
    let mut store = ChunkStore::new(9);
    store.get_or_generate(ChunkCoord::new(0, 0, 0), &NoiseGen).columns_mut();

    let clock = SimClock::restore(
        Tick(1_234_567),
        au_core::time::SimInstant { secs: 999, frac: 0xABCD },
        SimDuration::from_secs(3600),
    );
    let mut log = EventLog::new(10);
    log.emit(Event { tick: Tick(7), layer: Layer::Universe, kind: 1, a: 2, b: 3, c: 4 });

    let bytes = snapshot::encode(9, &clock, &store, &log, &[], [0; 8]);
    let snap = snapshot::decode(&bytes, &registry()).unwrap();

    assert_eq!(snap.clock.tick(), Tick(1_234_567));
    assert_eq!(snap.clock.now().frac, 0xABCD, "sub-second precision must survive a save");
    assert_eq!(snap.clock.scale(), SimDuration::from_secs(3600));
    assert_eq!(snap.events.len(), 1);
    assert_eq!(snap.events[0].a, 2);
}

#[test]
fn corruption_is_detected_not_absorbed() {
    let store = ChunkStore::new(1);
    let clock = SimClock::new(SimDuration::from_secs(1));
    let mut bytes = snapshot::encode(1, &clock, &store, &EventLog::new(10), &[], [0; 8]);

    let n = bytes.len();
    bytes[n / 2] ^= 0xFF;

    // A save file is a claim about the past. A damaged one must not load
    // "mostly".
    match snapshot::decode(&bytes, &registry()) {
        Err(CodecError::ChecksumMismatch { .. }) => {}
        Err(e) => panic!("wrong error: {}", e),
        Ok(_) => panic!("a corrupted save loaded anyway"),
    }
}

#[test]
fn generated_chunks_are_never_dirty() {
    // If a generator returned a dirty chunk, that chunk could never be evicted,
    // and the memory strategy would leak the entire world.
    let mut rng = Rng::derive(1, Domain::CHUNK_GEN, 0, 0);
    let c = NoiseGen.generate(ChunkCoord::new(0, 0, 0), &mut rng);
    assert!(!c.is_dirty());
    assert!(!EmptyGenerator.generate(ChunkCoord::new(0, 0, 0), &mut rng).is_dirty());
}
