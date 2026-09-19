//! Whole-universe invariants.
//!
//! You cannot unit-test emergence. There is no assertion of the form "a wolf
//! appears by tick 10^9" — there are no wolves, and there is no *should*. That
//! is the point of the project and it is also its central testing problem.
//!
//! What remains testable are **structural** properties, and they turn out to be
//! the ones that matter:
//!
//! * The same seed produces the same universe.
//! * A saved universe resumes into the same future.
//! * Nothing the engine does for its own convenience is visible to the world.
//!
//! If those three hold, then when Phase 8 produces no life after four hundred
//! million simulated years, we can *replay* it and find out why. If they don't,
//! we can only shrug. Everything downstream rests on this file.

use au_data::chunk::ChunkCoord;
use au_sim::{boot, Config, Simulation, DEFAULT_CONFIG};

fn cfg(seed: u64) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.set("world.seed", &seed.to_string());
    c
}

// ─── Determinism ─────────────────────────────────────────────────────────────

#[test]
fn the_same_seed_yields_the_same_universe() {
    let mut a = boot(Some(4242));
    let mut b = boot(Some(4242));
    a.run(200_000);
    b.run(200_000);
    assert_eq!(a.hash(), b.hash());
}

#[test]
fn different_seeds_yield_different_universes() {
    // Weak as a test, strong as a smoke alarm: if this ever fails, the seed is
    // not actually reaching the simulation and every "deterministic!" result
    // above is vacuous.
    let mut a = boot(Some(1));
    let mut b = boot(Some(2));
    a.run(1000);
    b.run(1000);
    assert_ne!(a.hash(), b.hash());
}

/// Determinism must hold *tick by tick*, not merely at the finish line.
///
/// Two runs that diverge and then coincidentally reconverge would pass an
/// end-state check. Comparing the whole trajectory is what lets us bisect to
/// the exact tick where a future bug enters.
#[test]
fn the_entire_history_matches_not_just_the_ending() {
    let trace = |seed: u64| {
        let mut s = boot(Some(seed));
        (0..5000).map(|_| { s.tick(); s.hash() }).collect::<Vec<_>>()
    };
    assert_eq!(trace(7), trace(7));
}

// ─── The engine must be invisible ────────────────────────────────────────────

/// The bug I actually wrote, now nailed down.
///
/// My first draft emitted a `CHUNK_EVICTED` event on eviction. Events feed the
/// world hash — so memory pressure would have become part of the identity of
/// the universe, and a machine with less RAM would have computed a *different
/// history*. It would have surfaced years from now as "evolution runs slightly
/// differently on the server".
///
/// Invariant: **nothing the engine does to manage memory may be observable to
/// the simulation.**
#[test]
fn memory_management_is_invisible_to_the_universe() {
    let mut hoarder = Simulation::new({
        let mut c = cfg(31337);
        c.set("chunks.evict_every_ticks", "0"); // never free anything
        c
    });
    let mut miser = Simulation::new({
        let mut c = cfg(31337);
        c.set("chunks.evict_every_ticks", "1"); // free everything, constantly
        c
    });

    for s in [&mut hoarder, &mut miser] {
        for x in -3..=3 {
            for z in -3..=3 {
                s.world.chunk(ChunkCoord::new(x, 0, z));
            }
        }
        s.run(100_000);
    }

    assert!(hoarder.world.chunks.stat_evicted == 0);
    assert!(miser.world.chunks.stat_evicted > 0, "the miser should have evicted");
    assert_eq!(
        hoarder.hash(),
        miser.hash(),
        "two hosts with different memory computed different universes"
    );
}

/// Likewise: *which parts of the world someone looked at* must not change what
/// the world is. Observation is not interaction.
#[test]
fn touching_chunks_does_not_change_history() {
    let mut untouched = boot(Some(88));
    let mut prodded = boot(Some(88));

    for i in 0..500 {
        prodded.world.chunk(ChunkCoord::new(i % 20, 0, i / 20));
    }
    untouched.run(50_000);
    prodded.run(50_000);

    assert_eq!(untouched.hash(), prodded.hash());
}

// ─── Persistence ─────────────────────────────────────────────────────────────

/// The test that says the past can be trusted.
///
/// A universe saved at tick K and resumed must produce, at tick N, exactly the
/// universe that never stopped. Any divergence means live state exists in RAM
/// that isn't in the snapshot — which means saved worlds are quietly *different*
/// worlds, and every family tree and lineage the player ever inspects is a
/// fabrication.
#[test]
fn a_saved_universe_resumes_into_the_same_future() {
    const N: u64 = 200_000;
    const K: u64 = 70_000;

    let mut straight = boot(Some(999));
    straight.run(N);

    let mut split = boot(Some(999));
    split.run(K);
    let bytes = split.save();
    let mut resumed = Simulation::load(&bytes, cfg(999)).unwrap();
    resumed.run(N - K);

    assert_eq!(resumed.world.clock.tick(), straight.world.clock.tick());
    assert_eq!(
        resumed.hash(),
        straight.hash(),
        "the resumed universe is not the one we saved"
    );
}

/// Deep time must survive a save. The clock crosses several rescale boundaries
/// here, so the *time scale itself* is part of what has to round-trip — if it
/// reset to 1 s/tick on load, the resumed world would age at the wrong rate and
/// nothing would ever notice until the ecologies disagreed.
#[test]
fn the_time_scale_survives_a_save() {
    let mut s = boot(Some(3));
    s.run(300_000); // several rescale boundaries at 65536 ticks
    let scale_before = s.world.clock.scale();
    assert!(scale_before.secs > 1, "deep-time acceleration should have engaged");

    let restored = Simulation::load(&s.save(), cfg(3)).unwrap();
    assert_eq!(restored.world.clock.scale(), scale_before);
    assert_eq!(restored.world.clock.now(), s.world.clock.now());
    assert_eq!(restored.hash(), s.hash());
}

// ─── Time ────────────────────────────────────────────────────────────────────

/// Deep-time acceleration (ARCHITECTURE §3): chemistry wants seconds, evolution
/// wants millennia, and the gap between them is ~10^14. No machine runs 10^14
/// ticks, so the clock coarsens as the universe ages.
#[test]
fn deep_time_acceleration_reaches_deep_time() {
    let mut s = boot(Some(1));
    s.run(2_000_000);

    let years = s.world.clock.now().as_years_f64();
    assert!(
        years > 1000.0,
        "2M ticks should have crossed millennia, not {:.2} years — \
         the accelerator is not engaging",
        years
    );

    // And it must stop, not run away to nonsense.
    let cap = s.world.config.u64_or("clock.rescale.max_seconds_per_tick", 0);
    assert!(s.world.clock.scale().secs <= cap, "the clock blew past its ceiling");

    // Every rescale is in the causal record. When something dies in a single
    // tick in Phase 9, this is the first thing to check.
    assert!(s.world.events.total() > 0, "rescales must be recorded as history");
}

#[test]
fn the_scheduler_is_honest_about_what_it_ran() {
    let mut s = boot(Some(1));
    s.run(65_536 * 4);
    let report = s.schedule_report();
    let (_, _, period, runs) = report.iter().find(|r| r.0 == "universe.time_scale").unwrap();
    assert_eq!(*period, 65_536);
    assert_eq!(*runs, 4, "a 65536-period system must run exactly 4 times in 4 periods");
}

// ─── The Golden Rule ─────────────────────────────────────────────────────────

/// Phase 1 ships an *empty* universe, and that is the correct state.
///
/// It would have been easy — and tempting, and fatal — to seed the world with
/// some placeholder terrain so the demo looked less bare. Matter does not exist
/// until Phase 3. Anything in the world right now would be content, and content
/// is the one thing the Golden Rule forbids. This test fails the day someone
/// gets impatient.
#[test]
fn the_world_contains_nothing_because_nothing_exists_yet() {
    let mut s = boot(Some(1));
    for x in -5..=5 {
        for z in -5..=5 {
            let c = s.world.chunk(ChunkCoord::new(x, 0, z));
            assert!(
                c.columns().is_empty(),
                "a chunk has content, but no layer that could produce content exists yet"
            );
        }
    }
    assert_eq!(s.world.entities.live_count(), 0, "there are no organisms in Phase 1");
}
