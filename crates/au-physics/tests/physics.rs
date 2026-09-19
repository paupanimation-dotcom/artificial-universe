//! Physics, checked against reality.
//!
//! # How do you test a universe that is supposed to surprise you?
//!
//! You cannot assert that anything *interesting* happens — that is the whole
//! point of the project. But physics is different from ecology: physics has
//! **right answers**, worked out over three centuries, and if this engine does
//! not reproduce them then every emergent thing built on top of it is emerging
//! from a lie.
//!
//! So these tests are of two kinds:
//!
//! * **Conservation** — bit-exact, non-negotiable, and the reason the state is
//!   integers. Evolution is an adversarial optimiser; a leak here becomes an
//!   organism in Phase 9.
//! * **Validation** — the engine is given real measured substances and asked to
//!   reproduce results that were known before computers existed: a linear
//!   conduction profile, a radiative equilibrium temperature, a latent-heat
//!   plateau of exactly the right duration. PROJECT_VISION allows precisely this
//!   ("Earth is only a scientific reference for how these systems work"). The
//!   substances are a *measuring instrument*, not content.

use au_physics::*;

const RHO_SILICATE: f64 = 3000.0;
const RHO_WATER: f64 = 1000.0;

/// Measured properties of real rock. A ruler, not a rock.
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
        // Rock melts at a *higher* temperature under pressure. Deep rock can
        // therefore be solid while shallower, cooler rock is molten — structure
        // that nobody put there.
        dtm_dp: 1.0e-7,
        dtb_dp: 1.0e-7,
        emissivity: 0.9,
        mu_liquid: 100.0,
        mu_gas: 1.0e-5,
        alpha: 3.0e-5,
    }
}

/// Water, with its genuinely strange sign. Ice melts at a *lower* temperature
/// under pressure — the one common substance that does. It is why ice floats,
/// why lakes freeze from the top down, and why anything survived the winter.
fn water() -> Material {
    Material {
        name: "water".into(),
        c_solid: 2093.0,
        c_liquid: 4184.0,
        c_gas: 1996.0,
        k_solid: 2.2,
        k_liquid: 0.6,
        k_gas: 0.025,
        melt_k: 273.15,
        boil_k: 373.15,
        latent_fusion: 333_550.0,
        latent_vapor: 2_257_000.0,
        dtm_dp: -7.4e-8,
        dtb_dp: 2.76e-4,
        emissivity: 0.96,
        mu_liquid: 1.0e-3,
        mu_gas: 1.0e-5,
        alpha: 2.1e-4,
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
    /// Fill a grid with one substance at one temperature. An *initial condition*.
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

// ═══ CONSERVATION — the non-negotiable one ═══════════════════════════════════

/// **The most important test in the project.**
///
/// Not "conserved to within 1e-9". Not "the error is small". *Exactly* conserved,
/// as an integer identity, over a hundred thousand substeps with heat pouring in
/// at the floor and radiating away at the sky:
///
/// ```text
/// E(now) − E(start)  ==  in − out
/// ```
///
/// If this ever fails by a single microjoule, stop everything. A leak is not a
/// rounding error, it is a **free-energy exploit**, and by Phase 9 something will
/// be alive that eats it.
#[test]
fn energy_is_conserved_exactly() {
    let spec = GridSpec::new(4, 24, 4, 10.0);
    let mut slab = Slab::uniform(spec, silicate(), RHO_SILICATE, 300.0);
    let b = Boundary {
        bottom_flux: 12.0,
        top_radiates: true,
        space_k: 2.7,
        surface_pressure: 0.0,
        gravity: 9.81,
        ..Default::default()
    };

    let e0 = slab.total();
    let mut ledger = Ledger::default();
    let mut scratch = Scratch::new();
    let reg = std::mem::take(&mut slab.reg);
    let mut v = slab.view();

    let mut total_substeps = 0u64;
    for _ in 0..2000 {
        let r = step(&mut v, &reg, &b, 5.0e6, 64, 0.4, &mut scratch);
        ledger.energy_in += r.energy_in;
        ledger.energy_out += r.energy_out;
        total_substeps += r.substeps as u64;
    }
    let e1 = v.total_energy();
    drop(v);

    assert!(total_substeps > 1000, "the test should have done real work");
    assert_eq!(
        e1 - e0,
        ledger.net_energy().0,
        "ENERGY LEAK. Off by {} µJ over {} substeps. This is not a rounding \
         error — it is a free-energy exploit waiting for something to evolve \
         into it.",
        (e1 - e0) - ledger.net_energy().0,
        total_substeps
    );
}

/// A closed box changes by *nothing*. Not nearly nothing. Nothing.
#[test]
fn a_closed_system_never_changes_its_total() {
    let spec = GridSpec::new(6, 6, 6, 50.0);
    let mut slab = Slab::uniform(spec, silicate(), RHO_SILICATE, 500.0);
    // A hot spot, so there is something to equilibrate and every face is busy.
    let hot = spec.idx(1, 1, 1);
    slab.energy[hot] *= 3;

    let e0 = slab.total();
    let b = Boundary::default(); // sealed: no source, no sky, no gravity
    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    let mut v = slab.view();
    for _ in 0..5000 {
        step(&mut v, &reg, &b, 1.0e6, 64, 0.4, &mut scratch);
    }
    let e1 = v.total_energy();
    assert_eq!(e0, e1, "a sealed box leaked {} µJ", e1 - e0);
}

/// And it must reach *equilibrium* — a single temperature — not merely conserve
/// its total while oscillating forever.
#[test]
fn a_closed_system_equilibrates() {
    let spec = GridSpec::new(8, 8, 1, 20.0);
    let mut slab = Slab::uniform(spec, silicate(), RHO_SILICATE, 400.0);
    slab.energy[spec.idx(0, 0, 0)] *= 4;

    let b = Boundary::default();
    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    {
        let mut v = slab.view();
        // One diffusion time for a 160 m slab of rock is L²/D ≈ 2.6e10 s.
        // Equilibration needs several of those. Physics does not care that we
        // are in a hurry.
        for _ in 0..12000 {
            step(&mut v, &reg, &b, 2.0e7, 64, 0.4, &mut scratch);
        }
    }
    slab.reg = reg;
    let t = slab.temps(&b);
    let (lo, hi) = t.iter().fold((f64::MAX, f64::MIN), |(l, h), &x| (l.min(x), h.max(x)));
    assert!(
        hi - lo < 1.0,
        "should have equilibrated; spread is still {:.2} K ({:.1}..{:.1})",
        hi - lo,
        lo,
        hi
    );
}

// ═══ VALIDATION — does it reproduce known physics? ═══════════════════════════

/// Fourier, 1822. At steady state, a constant heat flux through a slab produces
/// a **linear** temperature profile with slope q/k. If the engine gets this
/// wrong, it does not conduct heat; it does something that superficially
/// resembles conducting heat.
#[test]
fn steady_state_conduction_is_linear() {
    let spec = GridSpec::new(1, 32, 1, 10.0);
    let q = 0.1; // W/m²
    let mut slab = Slab::uniform(spec, silicate(), RHO_SILICATE, 40.0);
    let b = Boundary {
        bottom_flux: q,
        top_radiates: true,
        space_k: 2.7,
        surface_pressure: 0.0,
        gravity: 0.0,
        ..Default::default()
    };

    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    {
        let mut v = slab.view();
        for _ in 0..60000 {
            step(&mut v, &reg, &b, 2.0e7, 64, 0.4, &mut scratch);
        }
    }
    slab.reg = reg;
    let t = slab.temps(&b);

    // dT/dy = −q/k, so every cell is (q/k)·dx colder than the one below it.
    let expected_step = q / silicate().k_solid * spec.cell_m;
    for y in 1..spec.ny - 1 {
        let d = t[spec.idx(0, y - 1, 0)] - t[spec.idx(0, y, 0)];
        assert!(
            (d - expected_step).abs() < 0.02 * expected_step.abs().max(1e-3),
            "profile is not linear at y={}: step {:.4} K, Fourier says {:.4} K",
            y,
            d,
            expected_step
        );
    }
}

/// Stefan, 1879. A surface absorbing q W/m² and radiating freely settles at
///
/// ```text
/// T = (q / εσ)^(1/4)
/// ```
///
/// This equation is the reason planets have temperatures. It is, quite directly,
/// the reason the habitable zone is a zone. If the engine cannot reproduce it,
/// Phase 4 cannot have a climate.
#[test]
fn radiative_equilibrium_matches_stefan_boltzmann() {
    let spec = GridSpec::new(1, 8, 1, 5.0);
    let q = 200.0;
    let eps = silicate().emissivity;
    let mut slab = Slab::uniform(spec, silicate(), RHO_SILICATE, 300.0);
    let b = Boundary {
        bottom_flux: q,
        top_radiates: true,
        space_k: 2.7,
        surface_pressure: 0.0,
        gravity: 0.0,
        ..Default::default()
    };

    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    {
        let mut v = slab.view();
        for _ in 0..40000 {
            step(&mut v, &reg, &b, 1.0e6, 128, 0.4, &mut scratch);
        }
    }
    slab.reg = reg;
    let t = slab.temps(&b);
    let t_top = t[spec.idx(0, spec.ny - 1, 0)];

    let analytic = (q / (eps * SIGMA)).powf(0.25);
    assert!(
        (t_top - analytic).abs() / analytic < 0.01,
        "surface settled at {:.2} K; Stefan–Boltzmann says {:.2} K",
        t_top,
        analytic
    );
}

/// The plateau.
///
/// Heat ice at constant power and its temperature climbs — then **stops**, dead,
/// while 333.55 kJ/kg goes into breaking the crystal instead of into motion. Only
/// when the last of it has melted does the temperature start moving again.
///
/// That flat stretch is a thermostat with no thermostat in it. It is why Earth
/// has a climate rather than a temperature swing, and it is the first genuinely
/// *emergent regulation* in this codebase: nobody wrote a stabiliser, it falls
/// out of the thermodynamics.
///
/// The test asserts the plateau lasts exactly m·L/P seconds. Not approximately.
#[test]
fn latent_heat_produces_a_plateau_of_exactly_the_right_length() {
    let spec = GridSpec::new(1, 1, 1, 1.0); // one cell: no conduction, just heating
    let mut slab = Slab::uniform(spec, water(), RHO_WATER, 250.0); // ice, below freezing
    let q = 5000.0; // W/m²
    let b = Boundary { bottom_flux: q, top_radiates: false, ..Default::default() };

    let w = water();
    let kg = RHO_WATER * spec.cell_volume();
    let power = q * spec.face_area(); // W

    // THE ENGINE WAS RIGHT AND THIS TEST WAS WRONG.
    //
    // I first compared against `w.melt_k` — 273.15 K, the textbook number. The
    // plateau appeared at 273.1575 K and the test found nothing.
    //
    // Because there is no atmosphere in this box. `melt_k` is quoted at one
    // standard atmosphere, and this cell sits in a vacuum, and water's
    // Clausius–Clapeyron slope is *negative* — so removing 101 kPa of pressure
    // pushes its melting point 7.5 mK *up*. The engine had already applied that
    // and was reporting the melting point of ice in a vacuum, correctly, to four
    // decimal places, while the test insisted on the melting point of ice in a
    // room.
    //
    // A nice reminder that the useful failures are the ones where the simulation
    // knows more physics than the person testing it.
    let melt = w.melt_at(0.0);
    let t_to_melt = kg * w.c_solid * (melt - 250.0) / power;
    let t_melting = kg * w.latent_fusion / power; // the plateau, analytically

    let dt = 1.0;
    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    let mut trace: Vec<(f64, f64)> = Vec::new();
    {
        let mut v = slab.view();
        for i in 0..((t_to_melt + t_melting) * 1.6) as usize {
            step(&mut v, &reg, &b, dt, 64, 0.4, &mut scratch);
            let s = transport::states(&v, &reg, &b);
            trace.push((i as f64 * dt, s[0].temperature));
        }
    }
    slab.reg = reg;

    // Find the flat stretch: consecutive samples within a millikelvin of 273.15.
    let flat: Vec<f64> = trace
        .iter()
        .filter(|(_, t)| (t - melt).abs() < 1e-3)
        .map(|(s, _)| *s)
        .collect();
    assert!(!flat.is_empty(), "the temperature never paused — there is no latent heat");

    let observed = flat.last().unwrap() - flat.first().unwrap();
    assert!(
        (observed - t_melting).abs() / t_melting < 0.02,
        "plateau lasted {:.0} s; m·L/P says {:.0} s",
        observed,
        t_melting
    );

    // And it must climb again afterwards, as liquid — not stay stuck.
    let end = trace.last().unwrap().1;
    assert!(end > melt + 1.0, "never resumed warming: ended at {:.2} K", end);
}

/// Clausius and Clapeyron. Pressure moves the melting point — *up* for rock,
/// which is why Earth's inner core is solid despite being hotter than the liquid
/// outer core, and *down* for water, which is the anomaly that lets ice float and
/// lakes freeze from the top.
///
/// Structure from pressure. Nobody designed a core.
#[test]
fn pressure_moves_the_melting_point_and_the_sign_depends_on_the_substance() {
    let deep = 1.0e9; // 1 GPa

    let rock = silicate();
    assert!(
        rock.melt_at(deep) > rock.melt_at(P_REF) + 50.0,
        "rock must be harder to melt under pressure"
    );

    let ice = water();
    assert!(
        ice.melt_at(deep) < ice.melt_at(P_REF) - 50.0,
        "water is the anomaly: pressure must *lower* its melting point"
    );
}

/// Gravity, integrated. P = ρgh — the weight of everything above you.
#[test]
fn hydrostatic_pressure_is_the_weight_above() {
    let spec = GridSpec::new(1, 40, 1, 25.0);
    let mut slab = Slab::uniform(spec, silicate(), RHO_SILICATE, 300.0);
    let b = Boundary { gravity: 9.81, surface_pressure: 0.0, ..Default::default() };

    let reg = std::mem::take(&mut slab.reg);
    let v = slab.view();
    let p = transport::pressure_profile(&v, &b);
    drop(v);
    slab.reg = reg;

    // Bottom cell: everything above it, plus half of itself.
    let depth = (spec.ny as f64 - 0.5) * spec.cell_m;
    let analytic = RHO_SILICATE * 9.81 * depth;
    let got = p[spec.idx(0, 0, 0)];
    assert!(
        (got - analytic).abs() / analytic < 1e-6,
        "pressure at the floor is {:.0} Pa; ρgh says {:.0} Pa",
        got,
        analytic
    );
    assert!(p[spec.idx(0, spec.ny - 1, 0)] < p[spec.idx(0, 0, 0)], "pressure must rise with depth");
}

/// The second law. Heat does not flow uphill, ever, and no arrangement of cells
/// may produce a spontaneous temperature *inversion* — because a cell that could
/// get hotter than its neighbours for free is a Maxwell's demon, and a Maxwell's
/// demon in Phase 2 is a perpetual-motion organism in Phase 9.
#[test]
fn heat_never_flows_from_cold_to_hot() {
    let spec = GridSpec::new(2, 1, 1, 1.0);
    let mut slab = Slab::uniform(spec, silicate(), RHO_SILICATE, 300.0);
    slab.energy[1] *= 2; // cell 1 is hotter

    let b = Boundary::default();
    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    let (e_cold_0, e_hot_0) = (slab.energy[0], slab.energy[1]);
    {
        let mut v = slab.view();
        for _ in 0..200 {
            step(&mut v, &reg, &b, 1.0, 64, 0.4, &mut scratch);
        }
    }
    slab.reg = reg;

    assert!(slab.energy[0] > e_cold_0, "the cold cell must warm");
    assert!(slab.energy[1] < e_hot_0, "the hot cell must cool");

    // And it must never overshoot into an inversion, no matter how long it runs.
    let reg = std::mem::take(&mut slab.reg);
    {
        let mut v = slab.view();
        for _ in 0..100_000 {
            step(&mut v, &reg, &b, 1.0, 64, 0.4, &mut scratch);
            assert!(
                v.energy[1] >= v.energy[0] - 2,
                "the cold cell overtook the hot one — the scheme is unstable"
            );
        }
    }
    slab.reg = reg;
}

// ═══ ROBUSTNESS ══════════════════════════════════════════════════════════════

/// Nothing conducts through nothing.
#[test]
fn vacuum_does_nothing() {
    let spec = GridSpec::new(3, 3, 1, 1.0);
    let mut slab = Slab::uniform(spec, silicate(), RHO_SILICATE, 500.0);
    for i in 0..spec.cells() {
        slab.material[i] = MaterialId::VACUUM.0;
        slab.mass[i] = 0;
        slab.energy[i] = 0;
    }
    let before = slab.energy.clone();
    let b = Boundary { bottom_flux: 1e6, top_radiates: true, ..Default::default() };
    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    {
        let mut v = slab.view();
        for _ in 0..100 {
            let r = step(&mut v, &reg, &b, 1e6, 64, 0.4, &mut scratch);
            assert_eq!(r.energy_in, Energy::ZERO, "you cannot heat a vacuum");
        }
    }
    assert_eq!(slab.energy, before);
}

/// The engine must refuse to lie.
///
/// Ask it to advance a million years in one step and it *cannot* — diffusion has
/// a speed limit and no amount of wanting makes it faster. The correct behaviour
/// is not to take the step and produce garbage, and not to silently take a
/// smaller one and pretend: it is to integrate as far as it honestly can, mark
/// the result `unresolved`, and hand the caller the shortfall.
///
/// This is the discovery that Phase 1's deep-time accelerator cannot survive
/// contact with physics. Better to find it here than in Phase 8.
#[test]
fn the_engine_admits_when_it_cannot_resolve_the_physics() {
    let spec = GridSpec::new(1, 16, 1, 10.0);
    let mut slab = Slab::uniform(spec, silicate(), RHO_SILICATE, 300.0);
    let b = Boundary { bottom_flux: 1.0, top_radiates: true, ..Default::default() };
    let reg = std::mem::take(&mut slab.reg);
    let mut scratch = Scratch::new();
    let mut v = slab.view();

    let ok = step(&mut v, &reg, &b, 1.0e6, 64, 0.4, &mut scratch);
    assert!(!ok.unresolved);
    assert_eq!(ok.dt_done, 1.0e6, "a resolvable step must advance the full dt");

    // One million years, in one tick.
    let too_far = step(&mut v, &reg, &b, 3.15e13, 64, 0.4, &mut scratch);
    assert!(too_far.unresolved, "the engine claimed it could fast-forward diffusion");
    assert!(too_far.dt_done < 3.15e13, "it must report the shortfall, not hide it");
    assert_eq!(too_far.substeps, 64, "it must spend its whole budget trying");

    // And crucially: it did not explode.
    let t = transport::states(&v, &reg, &b);
    for s in &t {
        assert!(s.temperature.is_finite() && s.temperature >= 0.0, "unstable: {:?}", s);
    }
}

/// Same inputs, same universe. Physics must not introduce a single bit of
/// nondeterminism — no iteration-order dependence, no accumulated float drift in
/// the state.
#[test]
fn physics_is_bit_for_bit_reproducible() {
    let run = || {
        let spec = GridSpec::new(4, 8, 4, 10.0);
        let mut slab = Slab::uniform(spec, silicate(), RHO_SILICATE, 600.0);
        slab.energy[spec.idx(2, 3, 1)] *= 2;
        let b = Boundary { bottom_flux: 3.0, top_radiates: true, gravity: 9.81, ..Default::default() };
        let reg = std::mem::take(&mut slab.reg);
        let mut scratch = Scratch::new();
        let mut v = slab.view();
        for _ in 0..3000 {
            step(&mut v, &reg, &b, 2.0e5, 64, 0.4, &mut scratch);
        }
        drop(v);
        slab.energy
    };
    assert_eq!(run(), run());
}
