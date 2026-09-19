//! The foundations, tested.
//!
//! These are not unit tests in the "did I typo the getter" sense. Each one
//! pins down a property that, if it ever silently broke, would make the whole
//! project unbuildable — and would break in a way that looks like "evolution is
//! weird" a hundred thousand lines from now.

use au_core::*;
use std::collections::HashSet;

// ─── Time ────────────────────────────────────────────────────────────────────

/// The one that justifies the whole fixed-point time system.
///
/// 0.1 seconds is not representable in binary. In f64, adding it a million
/// times gives you *not* 100,000 — it gives you 100,000.00000133288. That error
/// is a divergence, and divergence is death here. Fixed point cannot drift,
/// because there is nothing to round.
#[test]
fn time_does_not_drift() {
    let step = SimDuration::from_secs_f64(0.1);
    let mut clock = SimClock::new(step);
    for _ in 0..1_000_000 {
        clock.advance();
    }

    // Exact: a million additions equals one multiplication. Bit for bit.
    let expected = step.mul_u64(1_000_000);
    assert_eq!(clock.now().secs, expected.secs);
    assert_eq!(clock.now().frac, expected.frac);

    // And for contrast, what f64 accumulation would have given us:
    let mut naive = 0.0f64;
    for _ in 0..1_000_000 {
        naive += 0.1;
    }
    assert_ne!(naive, 100_000.0, "if this ever passes, f64 stopped being f64");
    assert_eq!(clock.now().secs, 100_000, "fixed point lands exactly");
}

/// Deep time has to actually reach deep time.
#[test]
fn time_survives_billions_of_years() {
    let year = SimDuration::from_years(1.0);
    let a_lot = year.mul_u64(4_500_000_000);
    let years = SimInstant::ORIGIN.advanced_by(a_lot).as_years_f64();
    assert!(
        (years - 4.5e9).abs() < 1.0,
        "4.5 billion years should round-trip; got {}",
        years
    );
}

#[test]
fn clock_scale_is_simulation_state() {
    let mut c = SimClock::new(SimDuration::from_secs(1));
    c.advance();
    c.set_scale(SimDuration::from_secs(1000));
    c.advance();
    // 1 second, then 1000 seconds. The scale change is history, not a UI knob.
    assert_eq!(c.now().secs, 1001);
    assert_eq!(c.tick(), Tick(2));
}

// ─── Randomness ──────────────────────────────────────────────────────────────

#[test]
fn rng_is_a_pure_function_of_its_inputs() {
    let mut a = Rng::derive(42, Domain::MUTATION, 7, 99);
    let mut b = Rng::derive(42, Domain::MUTATION, 7, 99);
    for _ in 0..1000 {
        assert_eq!(a.next_u64(), b.next_u64());
    }
}

/// The property that makes parallelism safe: neighbouring keys must not produce
/// correlated streams. If chunk (0,0,0) and chunk (0,0,1) shared randomness,
/// the world would grow visible seams — and worse, mutations in adjacent
/// populations would correlate, quietly breaking evolution.
#[test]
fn rng_streams_are_independent() {
    let mut firsts = HashSet::new();
    for key in 0..10_000u64 {
        let mut r = Rng::derive(1, Domain::CHUNK_GEN, key, 0);
        firsts.insert(r.next_u64());
    }
    assert_eq!(firsts.len(), 10_000, "adjacent keys collided — streams are correlated");

    // Same key, different domain ⇒ different stream. Otherwise mutation and
    // behaviour would draw the same numbers.
    let mut m = Rng::derive(1, Domain::MUTATION, 5, 5);
    let mut b = Rng::derive(1, Domain::BEHAVIOR, 5, 5);
    assert_ne!(m.next_u64(), b.next_u64());
}

#[test]
fn rng_is_uniform_enough_to_trust() {
    let mut r = Rng::derive(9, Domain::TEST, 0, 0);
    let n = 200_000;
    let mut buckets = [0u32; 10];
    for _ in 0..n {
        let f = r.next_f64();
        assert!((0.0..1.0).contains(&f));
        buckets[(f * 10.0) as usize] += 1;
    }
    // Expect 20,000 per bucket. Allow 5%: we need "not obviously broken", not
    // a statistics paper.
    for (i, &b) in buckets.iter().enumerate() {
        assert!(
            (b as i64 - 20_000).abs() < 1_000,
            "bucket {} had {} samples; distribution is skewed",
            i,
            b
        );
    }

    // `below` must be unbiased — mutation picks alleles with it.
    let mut counts = [0u32; 3];
    for _ in 0..90_000 {
        counts[r.below(3) as usize] += 1;
    }
    for c in counts {
        assert!((c as i64 - 30_000).abs() < 800, "below(3) is biased: {:?}", counts);
    }
}

// ─── Layers ──────────────────────────────────────────────────────────────────

/// The Golden Rule, as an assertion.
///
/// ARCHITECTURE says a higher layer must never produce what a lower layer
/// should. In most projects that sentence lives in a document and quietly rots.
/// Here the engine refuses to start.
#[test]
fn a_system_may_not_read_upward() {
    let mut s: Scheduler<u32> = Scheduler::new();

    // Chemistry reading Physics: fine. Causes flow upward.
    assert!(s
        .register(
            SystemDesc::new("chem", Layer::Chemistry).reads(&[Layer::Physics, Layer::Matter]),
            |_: &mut u32| {},
        )
        .is_ok());

    // Physics reading Ecology: forbidden. This is the engine asking "what
    // content does the world need?" instead of "what happens next?".
    let err = s
        .register(
            SystemDesc::new("gravity_but_it_knows_about_wolves", Layer::Physics)
                .reads(&[Layer::Ecology]),
            |_: &mut u32| {},
        )
        .unwrap_err();

    assert!(matches!(err, ScheduleError::ReadsAbove { .. }));
    let msg = err.to_string();
    assert!(msg.contains("LAYER VIOLATION"), "{}", msg);
}

#[test]
fn a_system_may_read_its_own_layer() {
    let mut s: Scheduler<u32> = Scheduler::new();
    // Organisms see other organisms. That is not a violation — that is ecology
    // waiting to happen.
    assert!(s
        .register(
            SystemDesc::new("predation", Layer::Organisms).reads(&[Layer::Organisms]),
            |_: &mut u32| {},
        )
        .is_ok());
}

// ─── Scheduling ──────────────────────────────────────────────────────────────

/// ARCHITECTURE §3: different systems tick at wildly different rates. This is
/// the mechanism that makes deep time affordable.
#[test]
fn multi_rate_dispatch() {
    use std::cell::RefCell;
    use std::rc::Rc;

    let log: Rc<RefCell<Vec<&'static str>>> = Rc::new(RefCell::new(Vec::new()));
    let mut s: Scheduler<()> = Scheduler::new();

    for (name, period, layer) in [
        ("physics", 1u64, Layer::Physics),
        ("chemistry", 4, Layer::Chemistry),
        ("evolution", 16, Layer::Evolution),
    ] {
        let l = log.clone();
        s.must_register(
            SystemDesc::new(name, layer).every(period),
            move |_: &mut ()| l.borrow_mut().push(name),
        );
    }

    for t in 0..64u64 {
        s.run_tick(Tick(t), &mut ());
    }

    let l = log.borrow();
    assert_eq!(l.iter().filter(|&&n| n == "physics").count(), 64);
    assert_eq!(l.iter().filter(|&&n| n == "chemistry").count(), 16);
    assert_eq!(l.iter().filter(|&&n| n == "evolution").count(), 4);
}

/// Order must come from the layer stack, never from registration order or a
/// hash map. A world whose system order wobbles is a world that cannot be
/// replayed.
#[test]
fn execution_order_is_bottom_up_regardless_of_registration() {
    use std::cell::RefCell;
    use std::rc::Rc;

    let log: Rc<RefCell<Vec<Layer>>> = Rc::new(RefCell::new(Vec::new()));
    let mut s: Scheduler<()> = Scheduler::new();

    // Register in deliberately backwards order.
    for (name, layer) in [
        ("civ", Layer::Civilization),
        ("eco", Layer::Ecology),
        ("phys", Layer::Physics),
        ("chem", Layer::Chemistry),
    ] {
        let l = log.clone();
        s.must_register(SystemDesc::new(name, layer), move |_: &mut ()| {
            l.borrow_mut().push(layer)
        });
    }
    s.run_tick(Tick(0), &mut ());

    assert_eq!(
        *log.borrow(),
        vec![Layer::Physics, Layer::Chemistry, Layer::Ecology, Layer::Civilization],
        "systems must run bottom-up; causes before effects"
    );
}

#[test]
fn phase_staggers_equal_period_systems() {
    use std::cell::RefCell;
    use std::rc::Rc;
    let hits: Rc<RefCell<Vec<u64>>> = Rc::new(RefCell::new(Vec::new()));
    let mut s: Scheduler<u64> = Scheduler::new();
    let h = hits.clone();
    s.must_register(
        SystemDesc::new("staggered", Layer::Physics).every(10).phase(3),
        move |t: &mut u64| h.borrow_mut().push(*t),
    );
    for t in 0..30u64 {
        let mut cur = t;
        s.run_tick(Tick(t), &mut cur);
    }
    assert_eq!(*hits.borrow(), vec![3, 13, 23]);
}

// ─── Identity ────────────────────────────────────────────────────────────────

/// Family trees are a promised feature. A stale handle must be *detectably*
/// stale, not silently point at whoever moved into the dead organism's slot.
#[test]
fn recycled_slots_do_not_resurrect_the_dead() {
    let mut a = EntityAllocator::new();
    let ancestor = a.alloc();
    assert!(a.is_live(ancestor));

    a.free(ancestor);
    assert!(!a.is_live(ancestor));

    let newcomer = a.alloc();
    assert_eq!(newcomer.index, ancestor.index, "slot should be reused");
    assert_ne!(newcomer.generation, ancestor.generation, "but it is not the same being");
    assert!(!a.is_live(ancestor), "the dead stay dead");
    assert!(a.is_live(newcomer));
}
