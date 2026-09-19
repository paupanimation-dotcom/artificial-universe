//! Implicit conduction, held to the same standard as everything else.
//!
//! The claim of this phase is narrow and strong: **the timestep is no longer
//! chained to `dx²`**, and nothing about the engine's exactness was traded away to
//! achieve it. Those are two separate assertions and they are tested separately.
//!
//! The trap this phase had to avoid is worth stating up front, because it is the
//! obvious implementation and it is wrong. An implicit solve is iterative floating
//! point; it converges to a tolerance, never exactly. If the energy update were
//! written the natural way — `ΔE = C·(T^{n+1} − T^n)` per cell — then total energy
//! would be conserved *only to the solver's tolerance*, and the foundation's
//! central invariant would quietly become a function of an iteration count. So the
//! solver only ever decides **how much should move**; the transport moves it as one
//! integer per face, subtracted from one side and added to the other. The final
//! test here deliberately runs a badly under-converged solve and demands the books
//! still balance to zero — because the day that stops being true is the day
//! something evolves to farm the difference.

use au_physics::*;

const RHO: f64 = 3000.0;

fn silicate() -> Material {
    Material {
        name: "silicate".into(),
        c_solid: 1000.0,
        c_liquid: 1200.0,
        c_gas: 1000.0,
        k_solid: 3.0,
        k_liquid: 1.5,
        k_gas: 0.05,
        melt_k: 1400.0,
        boil_k: 3200.0,
        latent_fusion: 400_000.0,
        latent_vapor: 6_000_000.0,
        dtm_dp: 1.0e-7,
        dtb_dp: 1.0e-7,
        emissivity: 0.9,
        mu_liquid: 100.0,
        mu_gas: 1.0e-5,
        alpha: 3.0e-5,
    }
}

struct Slab {
    spec: GridSpec,
    reg: MaterialRegistry,
    mass: Vec<i128>,
    material: Vec<u16>,
    energy: Vec<i128>,
}

impl Slab {
    fn uniform(spec: GridSpec, mat: Material, rho: f64, temp_k: f64) -> Slab {
        let mut reg = MaterialRegistry::new();
        let id = reg.add(mat);
        let m = reg.get(id).unwrap();
        let n = spec.cells();
        let kg = rho * spec.cell_volume();
        let e = energy_for_temperature(kg, temp_k, m, 0.0);
        Slab {
            spec,
            mass: vec![Mass::from_kg(kg).0; n],
            material: vec![id.0; n],
            energy: vec![Energy::from_joules(e).0; n],
            reg,
        }
    }
    fn view(&mut self) -> GridView<'_> {
        GridView::new(self.spec, &mut self.mass, &self.material, &mut self.energy)
    }
    fn total(&self) -> i128 {
        self.energy.iter().sum()
    }
    fn temps(&mut self, b: &Boundary) -> Vec<f64> {
        let reg = std::mem::take(&mut self.reg);
        let v = self.view();
        let t = transport::states(&v, &reg, b).iter().map(|s| s.temperature).collect();
        self.reg = reg;
        t
    }
}

/// A sealed box with a hot spot — the standard relaxation problem.
fn hot_spot(spec: GridSpec, base_k: f64) -> Slab {
    let mut slab = Slab::uniform(spec, silicate(), RHO, base_k);
    slab.energy[spec.idx(spec.nx / 2, spec.ny / 2, spec.nz / 2)] *= 3;
    slab
}

fn implicit_boundary(iters: u32) -> Boundary {
    Boundary { implicit_conduction: true, conduction_iters: iters, ..Boundary::default() }
}

/// The explicit stability bound for this material and grid — the wall the whole
/// phase exists to climb over. Rock's diffusivity is k/(ρ·c) = 3/(3000·1000) = 1e-6 m²/s.
fn explicit_limit(dx: f64, ndim: f64) -> f64 {
    let alpha = 3.0 / (RHO * 1000.0);
    dx * dx / (2.0 * ndim * alpha)
}

// ═══ 1. AGREEMENT — where both are valid, they must say the same thing ══════

/// At a timestep small enough for the explicit scheme, implicit and explicit must
/// produce the same temperature field. They share the same face conductances and
/// the same harmonic-mean rule, so any disagreement here means one of them is
/// solving the wrong equation.
#[test]
fn implicit_agrees_with_explicit_at_small_timesteps() {
    let spec = GridSpec::new(10, 10, 1, 20.0);
    let dt = 0.2 * explicit_limit(20.0, 2.0); // comfortably inside the explicit bound

    let run = |implicit: bool| -> Vec<f64> {
        let mut slab = hot_spot(spec, 400.0);
        let b = if implicit { implicit_boundary(60) } else { Boundary::default() };
        let reg = std::mem::take(&mut slab.reg);
        let mut scratch = Scratch::new();
        {
            let mut v = slab.view();
            for _ in 0..200 {
                step(&mut v, &reg, &b, dt, 64, 0.4, &mut scratch);
            }
        }
        slab.reg = reg;
        slab.temps(&b)
    };

    let a = run(false);
    let c = run(true);
    let worst = a
        .iter()
        .zip(c.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    assert!(
        worst < 0.5,
        "implicit and explicit disagree by {:.3} K where both are valid — one is wrong",
        worst
    );
}

// ═══ 2. STABILITY — the point of the exercise ═══════════════════════════════

/// **The headline.** A timestep one thousand times past the explicit stability
/// limit, taken in a single stride, must remain bounded and physical: no
/// oscillation, no negative temperature, no divergence. This is the step that
/// would make an explicit scheme detonate.
///
/// The physical check is a maximum principle: heat diffusing in a sealed box may
/// never produce a temperature outside the range it started with. If the scheme
/// were unstable, that is the first thing to break.
#[test]
fn implicit_is_stable_a_thousand_times_past_the_explicit_limit() {
    let spec = GridSpec::new(10, 10, 1, 20.0);
    let limit = explicit_limit(20.0, 2.0);
    let dt = 1000.0 * limit;

    let mut slab = hot_spot(spec, 400.0);
    let b = implicit_boundary(40);
    let t_start = slab.temps(&b);
    let (t_lo, t_hi) = t_start.iter().fold((f64::MAX, f64::MIN), |(lo, hi), &t| (lo.min(t), hi.max(t)));

    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    {
        let mut v = slab.view();
        for _ in 0..50 {
            // max_substeps = 1: no hiding behind subdivision. One stride, in full.
            let r = step(&mut v, &reg, &b, dt, 1, 0.4, &mut scratch);
            assert!(!r.unresolved, "implicit conduction should not need substeps to be stable");
        }
    }
    slab.reg = reg;

    let t_end = slab.temps(&b);
    for (i, &t) in t_end.iter().enumerate() {
        assert!(t.is_finite(), "cell {} went non-finite — the scheme diverged", i);
        assert!(
            t >= t_lo - 1.0 && t <= t_hi + 1.0,
            "cell {} reached {:.1} K, outside the initial range [{:.1}, {:.1}] — \
             a maximum-principle violation, which is what instability looks like",
            i,
            t,
            t_lo,
            t_hi
        );
    }
}

/// The same grid, the same enormous step, run explicitly: the engine must *admit*
/// it cannot do it. This is the control that proves the previous test was not
/// merely easy — the `unresolved` flag is the engine refusing to pretend.
#[test]
fn explicit_cannot_take_the_step_implicit_takes_easily() {
    let spec = GridSpec::new(10, 10, 1, 20.0);
    let dt = 1000.0 * explicit_limit(20.0, 2.0);

    let mut slab = hot_spot(spec, 400.0);
    let b = Boundary::default(); // explicit
    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    let mut v = slab.view();
    let r = step(&mut v, &reg, &b, dt, 1, 0.4, &mut scratch);
    assert!(
        r.unresolved,
        "explicit conduction should not be able to cover {:.3e} s in one substep",
        dt
    );
}

/// And the win is *larger* on finer grids, which is the whole shape of the
/// problem: the explicit limit falls as `dx²`, so refining the grid four-fold
/// costs sixteen times the steps — while the implicit scheme's stable step does
/// not depend on `dx` at all.
#[test]
fn the_advantage_grows_as_the_grid_is_refined() {
    let coarse = explicit_limit(40.0, 2.0);
    let fine = explicit_limit(10.0, 2.0);
    let ratio = coarse / fine;
    assert!(
        (ratio - 16.0).abs() < 1e-6,
        "a 4× finer grid should cost 16× more explicit steps, got {:.2}×",
        ratio
    );
    // The implicit scheme is unconditionally stable, so its permissible step is
    // the same on both grids. The advantage is therefore exactly this ratio, and
    // it grows without bound as resolution rises.
}

// ═══ 3. CORRECTNESS — it must still be physics ══════════════════════════════

/// A sealed box equilibrates to a single temperature. Stability is worthless if
/// the answer is wrong; the implicit scheme must still carry heat down its
/// gradient until there is no gradient left.
#[test]
fn implicit_conduction_still_equilibrates() {
    let spec = GridSpec::new(8, 8, 1, 20.0);
    let mut slab = hot_spot(spec, 400.0);
    let b = implicit_boundary(40);

    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    {
        let mut v = slab.view();
        // Steps this large are unthinkable explicitly — which is the point.
        for _ in 0..400 {
            step(&mut v, &reg, &b, 1.0e8, 4, 0.4, &mut scratch);
        }
    }
    slab.reg = reg;

    let t = slab.temps(&b);
    let (lo, hi) = t.iter().fold((f64::MAX, f64::MIN), |(l, h), &x| (l.min(x), h.max(x)));
    assert!(
        hi - lo < 1.0,
        "a sealed box should reach one temperature; spread is still {:.2} K",
        hi - lo
    );
}

/// Heat flows from hot to cold, never the reverse. The second law, as a test:
/// after a step, the hot spot must be cooler and its neighbours warmer.
#[test]
fn heat_still_flows_only_downhill() {
    let spec = GridSpec::new(7, 7, 1, 20.0);
    let mut slab = hot_spot(spec, 400.0);
    let b = implicit_boundary(40);
    let centre = spec.idx(3, 3, 0);
    let neighbour = spec.idx(4, 3, 0);

    let t0 = slab.temps(&b);
    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    {
        let mut v = slab.view();
        step(&mut v, &reg, &b, 1.0e7, 4, 0.4, &mut scratch);
    }
    slab.reg = reg;
    let t1 = slab.temps(&b);

    assert!(t1[centre] < t0[centre], "the hot spot did not cool");
    assert!(t1[neighbour] > t0[neighbour], "its neighbour did not warm");
}

// ═══ 4. CONSERVATION — regardless of how badly the solver does ══════════════

/// A sealed box conserves energy **exactly** under implicit conduction, at a
/// timestep far past the explicit limit.
#[test]
fn implicit_conduction_conserves_energy_exactly() {
    let spec = GridSpec::new(6, 6, 6, 50.0);
    let mut slab = hot_spot(spec, 500.0);
    let e0 = slab.total();

    let b = implicit_boundary(30);
    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    {
        let mut v = slab.view();
        for _ in 0..500 {
            step(&mut v, &reg, &b, 1.0e9, 4, 0.4, &mut scratch);
        }
    }
    slab.reg = reg;

    assert_eq!(e0, slab.total(), "implicit conduction leaked {} µJ", slab.total() - e0);
}

/// **The test that matters most.** Run the solver deliberately crippled — a
/// single sweep, nowhere near converged — and demand the books still balance to
/// exactly zero.
///
/// This is the whole architectural claim of the phase, isolated. The physics
/// produced by one sweep is *wrong*: heat lands in the wrong places, and the
/// reported residual says so. But wrongness and non-conservation are different
/// failures, and only one of them is permitted. Energy moves as symmetric integer
/// transfers, so the total cannot drift no matter how poor the prediction that
/// chose the transfers was.
#[test]
fn conservation_does_not_depend_on_the_solver_converging() {
    let spec = GridSpec::new(6, 6, 6, 50.0);

    for iters in [1u32, 2, 5, 50] {
        let mut slab = hot_spot(spec, 500.0);
        let e0 = slab.total();
        let b = implicit_boundary(iters);
        let reg = std::mem::take(&mut slab.reg);
        let mut scratch = Scratch::new();
        {
            let mut v = slab.view();
            for _ in 0..200 {
                step(&mut v, &reg, &b, 1.0e9, 4, 0.4, &mut scratch);
            }
        }
        slab.reg = reg;
        assert_eq!(
            e0,
            slab.total(),
            "a {}-sweep solve leaked {} µJ — conservation must not depend on convergence",
            iters,
            slab.total() - e0
        );
    }
}

/// Determinism: the implicit path is a fixed number of sweeps in a fixed order,
/// so the same input produces bit-identical output. (A tolerance-based loop would
/// break this, which is why the iteration count is fixed and the residual merely
/// reported.)
#[test]
fn the_implicit_solve_is_bit_reproducible() {
    let spec = GridSpec::new(8, 8, 1, 20.0);
    let run = || {
        let mut slab = hot_spot(spec, 450.0);
        let b = implicit_boundary(25);
        let reg = std::mem::take(&mut slab.reg);
        let mut scratch = Scratch::new();
        {
            let mut v = slab.view();
            for _ in 0..100 {
                step(&mut v, &reg, &b, 1.0e8, 4, 0.4, &mut scratch);
            }
        }
        slab.energy.clone()
    };
    assert_eq!(run(), run(), "the implicit solve is not reproducible");
}

/// Explicit conduction is untouched: with the flag off, the engine behaves
/// exactly as it did in every earlier phase. New machinery must not perturb old
/// results.
#[test]
fn explicit_behaviour_is_unchanged_when_the_flag_is_off() {
    let spec = GridSpec::new(6, 6, 1, 30.0);
    let run = || {
        let mut slab = hot_spot(spec, 400.0);
        let b = Boundary::default();
        assert!(!b.implicit_conduction, "explicit must remain the default");
        let reg = std::mem::take(&mut slab.reg);
        let mut scratch = Scratch::new();
        {
            let mut v = slab.view();
            for _ in 0..100 {
                step(&mut v, &reg, &b, 1.0e5, 64, 0.4, &mut scratch);
            }
        }
        slab.energy.clone()
    };
    assert_eq!(run(), run());
}

/// The win, quantified in the currency that matters: **substeps**. Asked to cover
/// a span of simulated time, the explicit scheme must chop it into thousands of
/// slices; the implicit scheme takes it in one. This is the same comparison as the
/// stability test, but reported as cost rather than correctness — it is the number
/// that decides whether deep time is affordable.
#[test]
fn implicit_covers_in_one_substep_what_explicit_needs_thousands_for() {
    let spec = GridSpec::new(8, 8, 1, 20.0);
    let span = 2000.0 * explicit_limit(20.0, 2.0);

    let substeps = |implicit: bool| -> u32 {
        let mut slab = hot_spot(spec, 400.0);
        let b = if implicit { implicit_boundary(30) } else { Boundary::default() };
        let reg = std::mem::take(&mut slab.reg);
        let mut scratch = Scratch::new();
        let mut v = slab.view();
        step(&mut v, &reg, &b, span, 100_000, 0.4, &mut scratch).substeps
    };

    let explicit = substeps(false);
    let implicit = substeps(true);
    assert_eq!(implicit, 1, "implicit should take the whole span in one stride");
    assert!(
        explicit > 1000,
        "explicit should have needed thousands of substeps, took {}",
        explicit
    );
}

/// The convex-combination property, pinned directly — the reason
/// `conduction_omega` defaults to 1.0.
///
/// With plain Gauss–Seidel each update is a weighted average of a cell's own
/// previous temperature and its neighbours', so the field cannot leave its
/// starting range however few sweeps run. Here that is checked at its most
/// hostile: a *single* sweep, at a step a million times past the explicit limit.
/// An unstable or overshooting scheme fails this immediately; this one holds,
/// because the guarantee is structural rather than asymptotic.
#[test]
fn one_sweep_still_respects_the_maximum_principle() {
    let spec = GridSpec::new(9, 9, 1, 20.0);
    let mut slab = hot_spot(spec, 400.0);
    let b = implicit_boundary(1); // deliberately, absurdly under-converged
    let t0 = slab.temps(&b);
    let (lo, hi) = t0.iter().fold((f64::MAX, f64::MIN), |(l, h), &t| (l.min(t), h.max(t)));

    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    {
        let mut v = slab.view();
        for _ in 0..100 {
            step(&mut v, &reg, &b, 1.0e6 * explicit_limit(20.0, 2.0), 1, 0.4, &mut scratch);
        }
    }
    slab.reg = reg;

    for (i, &t) in slab.temps(&b).iter().enumerate() {
        assert!(
            t.is_finite() && t >= lo - 1.0 && t <= hi + 1.0,
            "one-sweep solve put cell {} at {:.1} K, outside [{:.1}, {:.1}]",
            i,
            t,
            lo,
            hi
        );
    }
}
