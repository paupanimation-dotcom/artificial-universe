//! Transport: how energy moves, and why it can never be lost on the way.
//!
//! # The two rules this file exists to enforce
//!
//! **1. Flux form.** Energy does not move by "each cell updating itself". It
//! moves across *faces*. Each face is visited exactly once, a single integer is
//! computed, and that same integer is subtracted from one side and added to the
//! other. Conservation is then not a property we test for and hope to keep — it
//! is an identity, true by construction, for the same reason that moving a coin
//! from one pocket to another cannot change how much money you have.
//!
//! **2. Gather, then scatter.** All fluxes are computed against the state at the
//! *start* of the substep, accumulated into a delta buffer, and only then
//! applied. Updating in place would mean cell (1,0) sees a neighbour that cell
//! (0,0) has already modified — making the result depend on iteration order. That
//! is a determinism bug, and it is the kind that survives for two years and then
//! surfaces as "the world is slightly different when we multithread it".
//!
//! # The uncomfortable discovery
//!
//! Diffusion has a speed limit. An explicit scheme is stable only while
//!
//! ```text
//! dt  ≤  dx² / (2 · ndim · D)
//! ```
//!
//! and if you exceed it the temperatures do not merely become inaccurate, they
//! oscillate to infinity within a few dozen steps.
//!
//! This collides head-on with the deep-time accelerator built in Phase 1. **You
//! cannot fast-forward a diffusion equation.** Doubling the seconds-per-tick does
//! not make heat spread twice as fast; it makes the arithmetic diverge. The only
//! honest options are to take more small steps (which costs exactly as much as
//! not fast-forwarding at all) or to *change the model*.
//!
//! So this file does the honest thing: it computes the largest timestep it can
//! actually integrate, sub-steps to stay inside it, and — when the clock demands
//! more than it can deliver — **caps the clock and says so**. The engine reports
//! that it cannot resolve the physics rather than quietly producing a number.
//!
//! That is a real constraint on the whole project and it is better to meet it in
//! Phase 2 than in Phase 8. The way out is not a faster integrator; it is level
//! of detail — deep time needs *coarse-grained physics*, not the same physics run
//! recklessly. Which is what LOD was always for.

use crate::fluid::{self, FluidScratch, FluidState};
use crate::grid::{Boundary, GridView};
use crate::material::{MaterialId, MaterialRegistry, SIGMA};
use crate::quantity::{Energy, Mass};
use crate::state::{derive, CellState};

/// Reusable working memory.
///
/// Owned by the caller and reused across ticks, because allocating four vectors
/// per chunk per tick would dominate the cost of the actual physics.
///
/// Everything in here is **derived**. None of it is simulation state, none of it
/// is snapshotted, and clearing it between ticks must change nothing. If it ever
/// does, something has smuggled state into the scratch space and determinism is
/// already gone.
#[derive(Default)]
pub struct Scratch {
    pressure: Vec<f64>,
    state: Vec<CellState>,
    /// Energy change for this substep, µJ. Sums to exactly (in − out).
    delta: Vec<i128>,
    /// Mass change for this substep, µg. Sums to exactly zero — nothing crosses
    /// the boundary in Phase 2b, matter only moves around inside.
    dmass: Vec<i128>,
    pub fluid: FluidScratch,
    pub implicit: crate::implicit::ImplicitScratch,
}

impl Scratch {
    pub fn new() -> Self {
        Self::default()
    }

    fn resize(&mut self, n: usize) {
        self.pressure.clear();
        self.pressure.resize(n, 0.0);
        self.state.clear();
        self.state.resize(n, CellState::VACUUM);
        self.delta.clear();
        self.delta.resize(n, 0);
        self.dmass.clear();
        self.dmass.resize(n, 0);
    }
}

/// What one call to [`step`] actually did.
#[derive(Clone, Copy, Debug, Default)]
pub struct StepReport {
    /// Energy that entered through a boundary. Exact.
    pub energy_in: Energy,
    /// Energy that left through a boundary. Exact.
    pub energy_out: Energy,
    /// How many stable substeps this tick needed.
    pub substeps: u32,
    /// The largest timestep this grid could be integrated with, in seconds.
    pub stable_dt: f64,
    /// **The engine admitting it could not do what was asked.** True when the
    /// requested dt exceeded what the substep budget could cover — the physics
    /// has been advanced by less than the clock asked for, and the caller must
    /// deal with that rather than pretend.
    pub unresolved: bool,
    /// Simulated seconds actually integrated. Equals the requested dt unless
    /// `unresolved`.
    pub dt_done: f64,
    pub max_temperature: f64,
    pub min_temperature: f64,

    /// Net external impulse delivered to the fluid this call, per axis. Exact.
    pub impulse: [i128; 3],
    /// Fastest thing in the grid, m·s⁻¹.
    pub max_speed: f64,
    /// Kinetic energy, J. Watch this to see an instability grow or die.
    pub kinetic_energy: f64,
    /// How badly the pressure solve failed to converge.
    pub projection_residual: f64,
    /// The worst ∇·u left in the field. The honest measure of how incompressible
    /// the fluid actually is, as opposed to how incompressible it was asked to be.
    pub max_divergence: f64,
}

/// Pressure at every cell. Public because the tests, the debugger and (later)
/// the renderer all need to see inside, and because a derived quantity nobody
/// can inspect is a derived quantity nobody can check.
pub fn pressure_profile(grid: &GridView, b: &Boundary) -> Vec<f64> {
    let mut out = vec![0.0; grid.spec.cells()];
    hydrostatic(grid, b, &mut out);
    out
}

/// Everything derivable about every cell, right now.
pub fn states(grid: &GridView, reg: &MaterialRegistry, b: &Boundary) -> Vec<CellState> {
    let p = pressure_profile(grid, b);
    let mut out = vec![CellState::VACUUM; grid.spec.cells()];
    derive_all(grid, reg, &p, &mut out);
    out
}

/// Hydrostatic pressure, integrated downward from the surface.
///
/// P(z) = P_surface + Σ ρ g dz. Pressure is *derived* — it is the weight of
/// everything above you — and it feeds straight back into the melting point,
/// which is why a deep enough layer can be solid while a shallower, cooler one
/// is molten. Nobody designs that. It just happens.
fn hydrostatic(grid: &GridView, b: &Boundary, out: &mut [f64]) {
    let spec = grid.spec;
    let v = spec.cell_volume();
    let dz = spec.cell_m;
    for z in 0..spec.nz {
        for x in 0..spec.nx {
            let mut p_above = b.surface_pressure;
            // Top down: y = ny-1 is the surface.
            for y in (0..spec.ny).rev() {
                let i = spec.idx(x, y, z);
                let rho = Mass(grid.mass[i]).as_kg() / v;
                let weight = rho * b.gravity * dz;
                // Pressure at the cell *centre* is the load above plus half its own.
                out[i] = p_above + 0.5 * weight;
                p_above += weight;
            }
        }
    }
}

fn derive_all(grid: &GridView, reg: &MaterialRegistry, pressure: &[f64], out: &mut [CellState]) {
    let v = grid.spec.cell_volume();
    for i in 0..grid.spec.cells() {
        let mid = grid.material_at(i);
        out[i] = match reg.get(mid) {
            Some(mat) if mid != MaterialId::VACUUM => derive(
                Mass(grid.mass[i]).as_kg(),
                Energy(grid.energy[i]).as_joules(),
                mat,
                pressure[i],
                v,
            ),
            _ => CellState::VACUUM,
        };
    }
}

/// The largest timestep this grid can be integrated with without diverging.
///
/// Two constraints, and we obey the tighter:
///
/// * **Diffusion:** dt ≤ dx² / (2 · ndim · D_max). Standard explicit-Euler
///   stability for the heat equation.
/// * **Radiation:** the boundary loss goes as T⁴, so its linearised response
///   time is C / (4εσA T³). A step longer than that will overshoot — the cell
///   radiates away more than it had, goes negative, and the whole grid rings.
///
/// Returns `f64::INFINITY` when nothing constrains it (an empty world imposes no
/// speed limit — which is exactly why Phase 1 could accelerate freely, and
/// exactly why it can no longer).
fn stable_dt(
    grid: &GridView,
    reg: &MaterialRegistry,
    states: &[CellState],
    b: &Boundary,
    safety: f64,
    fs: Option<&FluidState>,
) -> f64 {
    let spec = grid.spec;
    let ndim = spec.dimensionality().max(1);
    let dx = spec.cell_m;
    let area = spec.face_area();

    let mut d_max: f64 = 0.0;
    for s in states.iter() {
        if s.diffusivity > d_max {
            d_max = s.diffusivity;
        }
    }

    let mut dt = f64::INFINITY;
    // The dx² limit — and the whole reason implicit conduction exists. When the
    // conduction solve is implicit it is unconditionally stable, so this term
    // simply does not apply and the timestep is free of the grid spacing. What
    // remains below (radiation, advection, viscosity) is real and still binds.
    if d_max > 0.0 && !b.implicit_conduction {
        dt = dt.min(safety * dx * dx / (2.0 * ndim as f64 * d_max));
    }

    if b.top_radiates {
        for z in 0..spec.nz {
            for x in 0..spec.nx {
                let i = spec.idx(x, spec.ny - 1, z);
                let s = states[i];
                if s.heat_capacity <= 0.0 || s.temperature <= 0.0 {
                    continue;
                }
                let Some(mat) = reg.get(grid.material_at(i)) else { continue };
                let dloss_de = 4.0 * mat.emissivity * SIGMA * area * s.temperature.powi(3)
                    / s.heat_capacity;
                if dloss_de > 0.0 {
                    dt = dt.min(safety / dloss_de);
                }
            }
        }
    }

    // Fluids bring two more speed limits, and they stack on top of the thermal
    // one rather than replacing it.
    if let Some(fs) = fs {
        // Viscous: momentum diffuses, exactly like heat, and diverges past the
        // same bound. This is the one that makes real mantle convection
        // unreachable — the mantle's kinematic viscosity is ~10¹⁷ m²·s⁻¹, which
        // would demand dt ≈ 10⁻¹³ s. Explicit methods simply cannot go there.
        // Skipped entirely when the viscous term is solved implicitly: backward
        // Euler is unconditionally stable, so letting this bound survive would
        // impose the very limit the solver exists to remove — the same exemption
        // implicit conduction needed, and the same trap. (Phase 4b shipped a bug
        // where `stable_dt`'s infinity sentinel gained a second meaning and every
        // implicit run silently became a no-op; conservation tests passed
        // happily, because doing nothing conserves beautifully.)
        let nu_max = if b.viscous_iters > 0 {
            0.0
        } else {
            fs.nu.iter().cloned().fold(0.0f64, f64::max)
        };
        if nu_max > 0.0 {
            dt = dt.min(safety * dx * dx / (2.0 * ndim as f64 * nu_max));
        }
        // Advective (Courant): a parcel may not cross more than one cell in a
        // step, or the scheme is solving a different problem from the one it was
        // given.
        let u_max = fluid::max_speed(grid, b, fs);
        if u_max > 0.0 {
            dt = dt.min(safety * dx / u_max);
        }
    }
    dt
}

/// Advance the grid by `dt` simulated seconds.
///
/// Sub-steps internally to stay stable. If `max_substeps` is not enough, it
/// integrates as far as it honestly can, sets [`StepReport::unresolved`], and
/// reports the shortfall in `dt_done` — it does *not* take an unstable step and
/// it does *not* silently pretend the time passed.
pub fn step(
    grid: &mut GridView,
    reg: &MaterialRegistry,
    b: &Boundary,
    dt: f64,
    max_substeps: u32,
    safety: f64,
    scratch: &mut Scratch,
) -> StepReport {
    let n = grid.spec.cells();
    scratch.resize(n);

    let mut rep = StepReport { substeps: 0, dt_done: 0.0, ..Default::default() };
    if dt <= 0.0 || n == 0 {
        return rep;
    }

    let has_fluid = b.fluid && grid.fluid.is_some();

    hydrostatic(grid, b, &mut scratch.pressure);
    derive_all(grid, reg, &scratch.pressure, &mut scratch.state);
    let fs0 = if has_fluid {
        Some(fluid::fluid_state(grid, reg, &scratch.state))
    } else {
        None
    };
    let dt_stable = stable_dt(grid, reg, &scratch.state, b, safety, fs0.as_ref());
    rep.stable_dt = dt_stable;

    // An infinite stable step means *nothing constrains the timestep*. That used
    // to have exactly one cause — an inert grid with no diffusivity, nothing to
    // integrate — so returning early was a harmless shortcut.
    //
    // Implicit conduction gave the sentinel a second meaning: the stiff term is
    // now solved unconditionally, so it imposes no bound, and infinity can mean
    // "take the whole step in one stride" on a grid that very much does have work
    // to do. The old shortcut therefore skipped the entire physics step, silently
    // — every implicit run was a no-op, and only the physics tests caught it (the
    // conservation tests passed happily, because doing nothing conserves
    // beautifully).
    //
    // One stride is the right answer in *both* cases: an inert grid finds nothing
    // to do and costs nothing, and an implicit grid gets exactly the step it asked
    // for.
    let (substeps, h) = if !dt_stable.is_finite() {
        (1u32, dt)
    } else {
        let needed = (dt / dt_stable).ceil().max(1.0);
        let substeps = if needed > max_substeps as f64 {
            rep.unresolved = true;
            max_substeps
        } else {
            needed as u32
        };
        (substeps, if rep.unresolved { dt_stable } else { dt / substeps as f64 })
    };

    let spec = grid.spec;
    let area = spec.face_area();
    let dx = spec.cell_m;
    let mut e_in: i128 = 0;
    let mut e_out: i128 = 0;

    for s in 0..substeps {
        if s > 0 {
            derive_all(grid, reg, &scratch.pressure, &mut scratch.state);
        }

        // ── The fluid moves first, and only then is it allowed to carry things.
        //
        // Order matters and it is not arbitrary: momentum is advected using the
        // velocity field from the *start* of the substep, then projected to be
        // divergence-free, and only the projected field is used to advect mass
        // and energy. Advect scalars with a field that still has divergence in it
        // and matter accumulates in cells for no reason — you get spurious density
        // and the buoyancy it drives is an artefact.
        let fs = if has_fluid {
            let fs = fluid::fluid_state(grid, reg, &scratch.state);
            let imp = fluid::fluid_step(grid, b, &fs, h, &mut scratch.fluid);
            for a in 0..3 {
                rep.impulse[a] += imp[a];
            }
            Some(fs)
        } else {
            None
        };

        for d in scratch.delta.iter_mut() {
            *d = 0;
        }
        for d in scratch.dmass.iter_mut() {
            *d = 0;
        }

        // ── Conduction. Heat down its own gradient, whether or not anything flows.
        if b.implicit_conduction {
            // Backward Euler: unconditionally stable, so `h` may be anything the
            // caller asked for. The solve is approximate; the transport it emits
            // is exact. See `implicit.rs` for why those are different sentences.
            crate::implicit::conduct_implicit(
                spec,
                &scratch.state,
                area,
                dx,
                h,
                b.conduction_iters,
                b.conduction_omega,
                &mut scratch.implicit,
                &mut scratch.delta,
            );
        } else {
            conduct(spec, &scratch.state, area, dx, h, &mut scratch.delta);
        }

        // ── Advection. Heat, carried by the flow.
        if let Some(fs) = &fs {
            advect_energy(grid, b, fs, area, h, &mut scratch.delta);
        }

        // ── Boundaries. Every joule that crosses one is written down.
        let (bi, bo) = boundaries(grid, reg, &scratch.state, b, area, dx, h, &mut scratch.delta);
        e_in += bi;
        e_out += bo;

        for i in 0..n {
            grid.energy[i] += scratch.delta[i];
            grid.mass[i] += scratch.dmass[i]; // zero in Boussinesq; the channel Phase 3 will use
        }
        rep.substeps += 1;
    }

    rep.dt_done = h * substeps as f64;
    rep.energy_in = Energy(e_in);
    rep.energy_out = Energy(e_out);

    derive_all(grid, reg, &scratch.pressure, &mut scratch.state);
    rep.max_temperature = f64::MIN;
    rep.min_temperature = f64::MAX;
    for s in scratch.state.iter() {
        if s.heat_capacity > 0.0 {
            rep.max_temperature = rep.max_temperature.max(s.temperature);
            rep.min_temperature = rep.min_temperature.min(s.temperature);
        }
    }
    if rep.max_temperature == f64::MIN {
        rep.max_temperature = 0.0;
        rep.min_temperature = 0.0;
    }

    if has_fluid {
        let fs = fluid::fluid_state(grid, reg, &scratch.state);
        rep.max_speed = fluid::max_speed(grid, b, &fs);
        rep.kinetic_energy = fluid::kinetic_energy(grid, b, &fs);
        rep.projection_residual = scratch.fluid.residual();
        rep.max_divergence = scratch.fluid.max_divergence();
    }
    rep
}

/// Advection of energy by the (now divergence-free) velocity field.
///
/// Flux form and **upwind**. Upwind rather than central because central
/// differencing overshoots — it can hand a cell more energy than its neighbours
/// had between them, or take it below zero, and a cell with negative energy is a
/// cell with a negative temperature, from which there is no recovering. Upwind is
/// monotone: it cannot manufacture a new extremum. The price is numerical
/// diffusion of order u·dx/2, and it is worth paying.
///
/// # Why mass is *not* advected here, and what it cost me to learn that
///
/// I wrote this to advect mass as well, and the solver produced convection at
/// half the critical Rayleigh number — which is to say, it produced convection
/// out of nothing. With the perturbation set to *zero* it still stirred itself
/// into a 2 × 10⁻⁵ m·s⁻¹ circulation. That is not a fluid; that is a bug with
/// good manners.
///
/// Two things were wrong, and both are worth writing down.
///
/// **Boussinesq says ρ is constant.** That is not a detail of the approximation,
/// it *is* the approximation. Under a divergence-free velocity field with uniform
/// density, the net mass flux into every cell is exactly zero — so advecting mass
/// is a no-op in the continuum, and in the discrete it is a no-op plus rounding
/// error. I was integrating pure noise and calling it transport.
///
/// **And that noise was louder than the signal.** Temperature is derived as
/// E/(m·c). Jitter the mass by one quantum and the temperature moves with it. The
/// mass quantum is 1 µg in a 10⁻³ kg cell — a relative precision of 10⁻⁶ — while
/// the temperature difference driving the whole experiment was 1.9 × 10⁻⁴ K on a
/// 300 K base, a relative signal of 6.5 × 10⁻⁷.
///
/// **The signal was below the quantum.** I had asked the engine to resolve
/// something finer than the grain of its own state, and it answered with noise,
/// as it should have.
///
/// This is the price of exact integers, and it is a price worth naming: they buy
/// perfect conservation, and they charge for it in **dynamic range** — a hard,
/// knowable floor beneath which nothing is real. Floating point would not have
/// complained. It would have quietly leaked energy *and* quietly produced the
/// same fictitious convection, and I would have believed it.
///
/// Phase 3's chemical species will be advected by exactly this code path — a
/// concentration carried by a mass flux — and will inherit the same floor. Worth
/// knowing before, rather than after.
fn advect_energy(
    grid: &GridView,
    b: &Boundary,
    fs: &FluidState,
    area: f64,
    dt: f64,
    de: &mut [i128],
) {
    let t = fluid::Topo::new(grid.spec, b);
    let spec = grid.spec;
    let n = spec.cells();
    let uf = fluid::face_velocities(grid, b, fs);

    for a in 0..3 {
        if !t.active[a] {
            continue;
        }
        for i in 0..n {
            let c = t.coord(i);
            // The −a face of cell i, between cell (c − e_a) and cell c.
            let mut lo_c = c;
            lo_c[a] -= 1;
            let (Some(lo), Some(hi)) = (t.cell(lo_c), Some(i)) else { continue };
            if lo == hi {
                continue; // degenerate axis
            }
            let u = uf[a][i];
            if u == 0.0 {
                continue;
            }
            // Positive u flows from lo into hi.
            let up = if u > 0.0 { lo } else { hi };
            let m_up = Mass(grid.mass[up]).as_kg();
            if m_up <= 0.0 {
                continue;
            }
            let rho_up = m_up / spec.cell_volume();
            let e_specific = Energy(grid.energy[up]).as_joules() / m_up; // J/kg

            let mass_flux = rho_up * u * area * dt; // kg
            let energy_flux = mass_flux * e_specific; // J

            // One integer. Out of one cell, into the other. Conservation survives
            // the fluid intact.
            let qe = Energy::from_joules(energy_flux).0;
            de[lo] -= qe;
            de[hi] += qe;
        }
    }
}

/// Fourier's law, in flux form.
///
/// Face conductivity is the **harmonic** mean, not the arithmetic one. Two slabs
/// in series add their *resistances*, so an insulator next to a conductor
/// behaves like an insulator — which is the whole point of an insulator. The
/// arithmetic mean would let heat pour straight through a firewall, and the
/// resulting world would have no thermal structure at all.
fn conduct(
    spec: crate::grid::GridSpec,
    states: &[CellState],
    area: f64,
    dx: f64,
    dt: f64,
    delta: &mut [i128],
) {
    let mut face = |a: usize, bi: usize| {
        let (sa, sb) = (states[a], states[bi]);
        let (ka, kb) = (sa.conductivity, sb.conductivity);
        if ka <= 0.0 || kb <= 0.0 {
            return; // vacuum conducts nothing
        }
        let k = 2.0 * ka * kb / (ka + kb);
        let dtemp = sa.temperature - sb.temperature;
        if dtemp == 0.0 {
            return;
        }
        let watts = k * area * dtemp / dx;
        // One integer. Subtracted from one side, added to the other. This single
        // line is the entire conservation guarantee.
        let de = Energy::from_joules(watts * dt).0;
        delta[a] -= de;
        delta[bi] += de;
    };

    for z in 0..spec.nz {
        for y in 0..spec.ny {
            for x in 0..spec.nx {
                let i = spec.idx(x, y, z);
                if x + 1 < spec.nx {
                    face(i, spec.idx(x + 1, y, z));
                }
                if y + 1 < spec.ny {
                    face(i, spec.idx(x, y + 1, z));
                }
                if z + 1 < spec.nz {
                    face(i, spec.idx(x, y, z + 1));
                }
            }
        }
    }
}

/// Heat in at the floor; heat out at the sky. Returns (in, out) in µJ, exact.
#[allow(clippy::too_many_arguments)]
fn boundaries(
    grid: &GridView,
    reg: &MaterialRegistry,
    states: &[CellState],
    b: &Boundary,
    area: f64,
    dx: f64,
    dt: f64,
    delta: &mut [i128],
) -> (i128, i128) {
    let spec = grid.spec;
    let (mut e_in, mut e_out) = (0i128, 0i128);

    // ── Fixed-temperature walls. ────────────────────────────────────────────
    //
    // The condition that *defines* Rayleigh–Bénard: hold the floor hot, hold the
    // ceiling cold, and say nothing whatsoever about what the fluid should do
    // between them.
    //
    // Implemented as conduction to a ghost cell half a cell beyond the wall. The
    // energy this pushes in or pulls out is real energy crossing a real boundary,
    // so it is booked — a Dirichlet wall is an infinite reservoir, and an
    // unbooked infinite reservoir is a free lunch.
    for (row, temp) in [(0usize, b.bottom_temp), (spec.ny - 1, b.top_temp)] {
        let Some(t_wall) = temp else { continue };
        for z in 0..spec.nz {
            for x in 0..spec.nx {
                let i = spec.idx(x, row, z);
                let s = states[i];
                if s.heat_capacity <= 0.0 || s.conductivity <= 0.0 {
                    continue;
                }
                let watts = s.conductivity * area * (t_wall - s.temperature) / (0.5 * dx);
                let de = Energy::from_joules(watts * dt).0;
                delta[i] += de;
                if de >= 0 {
                    e_in += de;
                } else {
                    e_out += -de;
                }
            }
        }
    }

    // Floor: a fixed heat flux. Radiogenic heating, in a Phase 4 planet; a knob,
    // here.
    if b.bottom_flux != 0.0 {
        for z in 0..spec.nz {
            for x in 0..spec.nx {
                let i = spec.idx(x, 0, z);
                if states[i].heat_capacity <= 0.0 {
                    continue; // cannot heat a vacuum
                }
                let de = Energy::from_joules(b.bottom_flux * area * dt).0;
                delta[i] += de;
                // The books are kept with the *same integer* that was added to
                // the cell — not a recomputed float. Recomputing would introduce
                // exactly the kind of half-microjoule discrepancy that
                // conservation testing exists to catch.
                if de >= 0 {
                    e_in += de;
                } else {
                    e_out += -de;
                }
            }
        }
    }

    // Sky: Stefan–Boltzmann. This is the only reason a planet has a temperature
    // rather than a trend.
    if b.top_radiates {
        let sky4 = b.space_k.powi(4);
        for z in 0..spec.nz {
            for x in 0..spec.nx {
                let i = spec.idx(x, spec.ny - 1, z);
                let s = states[i];
                if s.heat_capacity <= 0.0 {
                    continue;
                }
                let Some(mat) = reg.get(grid.material_at(i)) else { continue };
                let watts = mat.emissivity * SIGMA * area * (s.temperature.powi(4) - sky4);
                let de = Energy::from_joules(watts * dt).0;
                delta[i] -= de;
                if de >= 0 {
                    e_out += de;
                } else {
                    // Colder than space. Absorbing from the sky. Physical, if odd.
                    e_in += -de;
                }
            }
        }
    }

    (e_in, e_out)
}
