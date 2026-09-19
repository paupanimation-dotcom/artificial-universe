//! The fluid solver.
//!
//! This is the hardest code in the project so far, and almost all of the
//! difficulty is in two places that look like implementation detail and are not.
//!
//! # 1. Incompressibility is a *global* constraint
//!
//! In an incompressible fluid, pressure is felt everywhere instantly. Push here
//! and the whole domain knows. That is not a local rule you can apply cell by
//! cell — it is an elliptic problem, and it must be *solved*, every step, over
//! the entire grid.
//!
//! There is no way around it. You can pretend the fluid is compressible and skip
//! the solve, but then the timestep is tied to the speed of sound — 5,770 m·s⁻¹
//! in rock, giving dt ≈ 2 ms — and the four billion years we eventually need are
//! not going to arrive.
//!
//! So: [`project`] solves a Poisson equation, every substep, with red-black SOR.
//! It costs more than everything else in this file combined.
//!
//! # 2. The checkerboard
//!
//! Put velocity and pressure both at cell centres — the obvious thing — and the
//! discrete pressure Laplacian only couples cells *two* apart. The odd and even
//! cells decouple entirely. The solver then converges happily to a pressure field
//! with an invisible checkerboard riding on it, which drives spurious flow that
//! looks, to the naked eye, exactly like convection.
//!
//! People have published convection cells that were checkerboard artefacts.
//!
//! Harlow & Welch's 1965 fix: velocities live on **faces**, pressure at centres.
//! The divergence of a cell then uses its own faces, the gradient at a face uses
//! its own cells, and the Laplacian is compact. It makes every loop in this file
//! fiddlier, and it is not optional.
//!
//! # What survives from Phase 2
//!
//! Everything. Advection and viscosity are still **symmetric flux** — one
//! integer, subtracted from one control volume and added to the other — so
//! momentum is conserved to the bit under internal exchange. Newton's third law
//! as an identity rather than an aspiration.
//!
//! What is *not* an internal exchange — buoyancy, the pressure force on the
//! walls, friction against a no-slip surface, momentum absorbed by a solid — gets
//! booked in the ledger. Which gives the identity this solver lives or dies by:
//!
//! ```text
//! Σ momentum(now) − Σ momentum(start)  ==  impulse     EXACTLY
//! ```
//!
//! Break it and the fluid can push against nothing, and something will
//! eventually evolve to swim on that.

use crate::grid::{Boundary, GridSpec, GridView, Wall};
use crate::material::{MaterialRegistry, Phase};
use crate::quantity::{Mass, Momentum};
use crate::state::CellState;

/// Grid topology: dimensions plus which axes wrap.
#[derive(Clone, Copy)]
pub struct Topo {
    pub spec: GridSpec,
    pub periodic: [bool; 3],
    /// An axis with one cell has no faces along it. Skip it entirely, or a cell
    /// becomes its own neighbour and the fluxes are nonsense.
    pub active: [bool; 3],
}

impl Topo {
    pub fn new(spec: GridSpec, b: &Boundary) -> Topo {
        let n = [spec.nx, spec.ny, spec.nz];
        Topo {
            spec,
            periodic: [
                b.x_wall.is_periodic(),
                b.y_wall.is_periodic(),
                b.z_wall.is_periodic(),
            ],
            active: [n[0] > 1, n[1] > 1, n[2] > 1],
        }
    }

    #[inline]
    pub fn n(&self, a: usize) -> isize {
        [self.spec.nx, self.spec.ny, self.spec.nz][a] as isize
    }

    /// Cell index, wrapping on periodic axes and returning `None` off a wall.
    #[inline]
    pub fn cell(&self, c: [isize; 3]) -> Option<usize> {
        let mut o = [0usize; 3];
        for a in 0..3 {
            let n = self.n(a);
            let mut v = c[a];
            if v < 0 || v >= n {
                if !self.periodic[a] || !self.active[a] {
                    return None;
                }
                v = v.rem_euclid(n);
            }
            o[a] = v as usize;
        }
        Some(self.spec.idx(o[0], o[1], o[2]))
    }

    #[inline]
    pub fn coord(&self, i: usize) -> [isize; 3] {
        let nx = self.spec.nx;
        let ny = self.spec.ny;
        [
            (i % nx) as isize,
            ((i / nx) % ny) as isize,
            (i / (nx * ny)) as isize,
        ]
    }
}

#[inline]
fn step_by(mut c: [isize; 3], a: usize, d: isize) -> [isize; 3] {
    c[a] += d;
    c
}

/// Everything the fluid step needs to know about a cell, derived fresh.
pub struct FluidState {
    /// Reference density, kg·m⁻³. Boussinesq: one number for the whole grid.
    ///
    /// This is the approximation, stated plainly. Density variations are assumed
    /// small and to matter *only* where they multiply gravity. It is what every
    /// mantle and ocean code on Earth does, and it is the assumption under which
    /// Ra_c = 657.511 is a true statement — so a solver that wants to be checked
    /// against that number has to make it.
    pub rho0: f64,
    /// Mean temperature. Buoyancy is measured against it, which makes the total
    /// buoyant force sum to exactly zero: a uniformly warm box does not levitate.
    pub t_ref: f64,
    /// Rigid cells. A solid does not flow — it is not a very thick liquid, and
    /// pretending otherwise would demand a viscous timestep of ~10⁻¹³ s.
    pub solid: Vec<bool>,
    /// Kinematic viscosity per cell, m²·s⁻¹.
    pub nu: Vec<f64>,
    /// Thermal expansion per cell, K⁻¹.
    pub alpha: Vec<f64>,
    pub temperature: Vec<f64>,
}

pub fn fluid_state(
    grid: &GridView,
    reg: &MaterialRegistry,
    states: &[CellState],
) -> FluidState {
    let n = grid.spec.cells();
    let v = grid.spec.cell_volume();
    let mut solid = vec![true; n];
    let mut nu = vec![0.0; n];
    let mut alpha = vec![0.0; n];
    let mut temperature = vec![0.0; n];

    let mut mass_sum = 0.0;
    let mut vol_sum = 0.0;
    let mut t_sum = 0.0;
    let mut fluid_cells = 0.0;

    for i in 0..n {
        let s = states[i];
        temperature[i] = s.temperature;
        let Some(mat) = reg.get(grid.material_at(i)) else { continue };
        if s.heat_capacity <= 0.0 {
            continue;
        }
        match mat.viscosity(s.phase) {
            None => {} // solid: stays rigid
            Some(mu) => {
                let rho = Mass(grid.mass[i]).as_kg() / v;
                if rho > 0.0 && mu > 0.0 {
                    solid[i] = false;
                    nu[i] = mu / rho;
                    alpha[i] = mat.alpha;
                    mass_sum += Mass(grid.mass[i]).as_kg();
                    vol_sum += v;
                    t_sum += s.temperature;
                    fluid_cells += 1.0;
                }
            }
        }
    }

    FluidState {
        rho0: if vol_sum > 0.0 { mass_sum / vol_sum } else { 0.0 },
        t_ref: if fluid_cells > 0.0 { t_sum / fluid_cells } else { 0.0 },
        solid,
        nu,
        alpha,
        temperature,
    }
}

/// Is the face on axis `a` at face-coordinate `c` a wall?
///
/// A wall is a face the fluid cannot cross: the edge of a walled domain, the
/// boundary of a solid, or the surface of a vacuum. **A solid is just a wall that
/// happens to be inside the grid** — which is how a convecting magma layer under
/// a frozen lid works without a single line of special-case code.
#[inline]
fn is_wall(t: &Topo, fs: &FluidState, a: usize, c: [isize; 3]) -> bool {
    if !t.active[a] {
        return true;
    }
    let hi = t.cell(c);
    let lo = t.cell(step_by(c, a, -1));
    match (lo, hi) {
        (Some(l), Some(h)) => fs.solid[l] || fs.solid[h],
        _ => true, // off the edge of a walled domain
    }
}

/// Velocity on the face at face-coordinate `c`, axis `a`. m·s⁻¹.
#[inline]
fn u_face(t: &Topo, fs: &FluidState, mom: &[&[i128]; 3], a: usize, c: [isize; 3]) -> f64 {
    if is_wall(t, fs, a, c) {
        return 0.0;
    }
    // The far face of the last cell isn't stored; on a periodic axis it *is*
    // face 0, which `Topo::cell` wraps to.
    match t.cell(c) {
        Some(i) => {
            let m = Momentum(mom[a][i]).as_si();
            let cell_mass = fs.rho0 * t.spec.cell_volume();
            if cell_mass > 0.0 {
                m / cell_mass
            } else {
                0.0
            }
        }
        None => 0.0,
    }
}

/// One fluid substep. Returns the external impulse delivered, per axis (µ-units).
#[allow(clippy::too_many_arguments)]
pub fn fluid_step(
    grid: &mut GridView,
    b: &Boundary,
    fs: &FluidState,
    dt: f64,
    scratch: &mut FluidScratch,
) -> [i128; 3] {
    let t = Topo::new(grid.spec, b);
    let n = grid.spec.cells();
    let dx = grid.spec.cell_m;
    let area = grid.spec.face_area();
    let vol = grid.spec.cell_volume();
    let cell_mass = fs.rho0 * vol;

    if fs.rho0 <= 0.0 {
        return [0; 3];
    }

    scratch.resize(n);
    let mut impulse = [0i128; 3];

    // ── Enforce the walls. ───────────────────────────────────────────────────
    // Any momentum sitting on a face that is now a wall — because the fluid
    // froze there, or because it was never legal — is absorbed. Absorbed, not
    // deleted: a wall that silently swallows momentum is a wall you can push off
    // for free.
    {
        let f = grid.fluid.as_mut().unwrap();
        let moms: [&mut [i128]; 3] = [f.mom_x, f.mom_y, f.mom_z];
        for (a, m) in moms.into_iter().enumerate() {
            for i in 0..n {
                let c = t.coord(i);
                if is_wall(&t, fs, a, c) && m[i] != 0 {
                    impulse[a] -= m[i];
                    m[i] = 0;
                }
            }
        }
    }

    // ── Gather the velocity field. ───────────────────────────────────────────
    let (mx, my, mz) = {
        let f = grid.fluid.as_ref().unwrap();
        (f.mom_x.to_vec(), f.mom_y.to_vec(), f.mom_z.to_vec())
    };
    let mom_ro: [&[i128]; 3] = [&mx, &my, &mz];

    for d in scratch.delta.iter_mut() {
        d.iter_mut().for_each(|v| *v = 0);
    }

    // ── Advection. Momentum carried by the flow. ─────────────────────────────
    //
    // The control volume for the a-momentum sits *on* the a-face, straddling two
    // cells. Its own faces are therefore offset by half a cell, and the velocity
    // advecting across them is the average of the two staggered velocities that
    // bracket it:
    //
    //     u_b at the +b face of CV(a, c)  =  ½[ u_b(c + e_b) + u_b(c + e_b − e_a) ]
    //
    // which — pleasingly — reduces to the cell-centred velocity when a == b, so
    // one expression covers both the along-axis and cross-axis cases.
    //
    // Upwind, not central. Central differencing is more accurate and can produce
    // negative energy, and a cell with negative energy is a cell with a negative
    // temperature, and there is no recovering from that.
    for a in 0..3 {
        if !t.active[a] {
            continue;
        }
        for i in 0..n {
            let c = t.coord(i);
            if is_wall(&t, fs, a, c) {
                continue;
            }
            let u_a_here = u_face(&t, fs, &mom_ro, a, c);

            for bx in 0..3 {
                if !t.active[bx] {
                    continue;
                }
                let cn = step_by(c, bx, 1); // the +b neighbouring CV
                if is_wall(&t, fs, a, cn) {
                    continue; // nothing on the other side to trade with
                }
                let Some(_) = t.cell(cn) else { continue };

                let u_b = 0.5
                    * (u_face(&t, fs, &mom_ro, bx, step_by(c, bx, 1))
                        + u_face(&t, fs, &mom_ro, bx, step_by(step_by(c, bx, 1), a, -1)));
                if u_b == 0.0 {
                    continue;
                }
                let u_a_there = u_face(&t, fs, &mom_ro, a, cn);
                let upwind = if u_b > 0.0 { u_a_here } else { u_a_there };

                // mass flux × the velocity it carries × dt
                let flux = fs.rho0 * u_b * area * upwind * dt;
                let q = Momentum::from_si(flux).0;

                // One integer. Out of one control volume, into the other.
                let ii = t.cell(c).unwrap();
                let jj = t.cell(cn).unwrap();
                scratch.delta[a][ii] -= q;
                scratch.delta[a][jj] += q;
            }
        }
    }

    // ── Viscosity. Momentum diffusing down its own gradient. ─────────────────
    for a in 0..3 {
        if !t.active[a] {
            continue;
        }
        // ── The implicit viscous solve, when it comes, goes here. ────────────
        //
        // This loop is the whole of explicit viscosity: for each face, a
        // symmetric exchange proportional to the velocity difference, with a
        // harmonic-mean viscosity because two slabs in series add resistances.
        // It is stable only while `dt <= dx²/(2·ndim·ν)`, and that bound is what
        // makes a micron-scale grid unaffordable — a 24x12 protocell convection
        // cell substeps for eleven minutes without finishing 460 ticks, which is
        // currently the blocker on protocell mobility and therefore on washout
        // and therefore on selection.
        //
        // **The change is far smaller than the roadmap assumed**, and the reason
        // is Phase 4b's discipline: *the solver proposes, the transport
        // disposes*. Implicit conduction did not rewrite its transport; it
        // computed a better field first and then moved exact integers across
        // faces exactly as before. The identical split works here:
        //
        //   1. solve `(I − dt·ν·∇²)·u_new = u_old` per component, Gauss–Seidel,
        //      `omega = 1.0` — plain, not over-relaxed, because the update is a
        //      convex combination and over-relaxation forfeits the maximum
        //      principle (Phase 4b found temperatures of 4244 K in a box that
        //      started between 400 K and 1200 K learning that);
        //   2. run **this loop unchanged**, reading `u_new` instead of `u_old`.
        //
        // Step 2 is a one-word change. Conservation therefore does not depend on
        // the solve converging: an under-converged field puts momentum in the
        // wrong place — honestly wrong, and visible in the residual — but cannot
        // create or destroy any, because the two halves of a symmetric integer
        // transfer are the same integer. Do **not** write the update as
        // `Δp = m·(u_new − u_old)` per cell; that is conservative only on a
        // perfectly converged solve, and it silently couples the engine's
        // central invariant to an iteration count.
        //
        // The one real difficulty is that velocity lives on faces (MAC), so each
        // component needs its own Laplacian assembly with its own wall handling
        // — `is_wall` and `nu_at_face` already answer both questions, and
        // `au-physics/src/implicit.rs` says at its head that it was written to
        // be reusable for exactly this.
        // ── Propose. ────────────────────────────────────────────────────────
        //
        // **Measured, and smaller than hoped.** On the 24x12 protocell
        // convection cell this buys 1.7x (7.8s -> 4.6s for 20 ticks), not the
        // orders of magnitude implicit conduction bought for heat. Viscosity was
        // therefore never the dominant cost there. What remains is the advective
        // Courant bound — `dt <= dx/u`, which no implicit viscous scheme touches
        // because advection is hyperbolic, not diffusive — and the pressure
        // projection's sixty sweeps per substep. Both are structural, and both
        // are worth measuring before anything else is built on the assumption
        // that fine-grid flow is now affordable. It is cheaper; it is not cheap.
        //
        // With `viscous_iters > 0`, solve the backward-Euler system first and
        // let the flux loop below read *that* field instead of the current one.
        // The loop is otherwise untouched, which is the whole point: transport
        // stays one exact integer per face, so an under-converged solve puts
        // momentum in the wrong place — honestly wrong — and cannot lose any.
        let implicit = b.viscous_iters > 0;
        if implicit {
            scratch.u0[a].clear();
            for i in 0..n {
                scratch.u0[a].push(u_face(&t, fs, &mom_ro, a, t.coord(i)));
            }
            let (u0, u) = (&scratch.u0[a], &mut scratch.u_new[a]);
            viscous_solve(&t, fs, a, dt, dx, b.viscous_iters, u0, u);
        }
        let u_at = |t: &Topo, fs: &FluidState, a: usize, c: [isize; 3]| -> f64 {
            if implicit {
                match t.cell(c) {
                    Some(i) if !is_wall(t, fs, a, c) => scratch.u_new[a][i],
                    _ => 0.0,
                }
            } else {
                u_face(t, fs, &mom_ro, a, c)
            }
        };

        for i in 0..n {
            let c = t.coord(i);
            if is_wall(&t, fs, a, c) {
                continue;
            }
            let ii = t.cell(c).unwrap();
            let u_here = u_at(&t, fs, a, c);
            let nu_here = nu_at_face(&t, fs, a, c);
            if nu_here <= 0.0 {
                continue;
            }

            for bx in 0..3 {
                if !t.active[bx] {
                    continue;
                }
                let cn = step_by(c, bx, 1);
                if !is_wall(&t, fs, a, cn) {
                    // Interior: a symmetric exchange, exactly conserving.
                    let Some(jj) = t.cell(cn) else { continue };
                    let u_there = u_at(&t, fs, a, cn);
                    let nu_there = nu_at_face(&t, fs, a, cn);
                    let nu = 2.0 * nu_here * nu_there / (nu_here + nu_there).max(1e-300);
                    let force = fs.rho0 * nu * area * (u_there - u_here) / dx;
                    let q = Momentum::from_si(force * dt).0;
                    scratch.delta[a][ii] += q;
                    scratch.delta[a][jj] -= q;
                } else if a != bx && b.wall(bx) == Wall::NoSlip {
                    // A no-slip wall drags on the flow along it. The wall keeps
                    // what it takes — book it, or the fluid gets free traction.
                    let force = fs.rho0 * nu_here * area * (0.0 - u_here) / (0.5 * dx);
                    let q = Momentum::from_si(force * dt).0;
                    scratch.delta[a][ii] += q;
                    impulse[a] += q;
                }
                // Free-slip: ∂u/∂n = 0 at the wall. No shear, nothing to book.

                // And the same on the −b side, so the CV is stressed on both faces.
                let cp = step_by(c, bx, -1);
                if is_wall(&t, fs, a, cp) && a != bx && b.wall(bx) == Wall::NoSlip {
                    let force = fs.rho0 * nu_here * area * (0.0 - u_here) / (0.5 * dx);
                    let q = Momentum::from_si(force * dt).0;
                    scratch.delta[a][ii] += q;
                    impulse[a] += q;
                }
            }
        }
    }

    // ── Buoyancy. The whole reason anything moves. ───────────────────────────
    //
    // ρ = ρ₀(1 − α(T − T̄)), so a parcel hotter than the mean is *lighter* than
    // the fluid around it, and the surrounding pressure gradient — which is
    // holding up fluid of density ρ₀ — over-supports it. Up it goes.
    //
    // Measured against the domain mean, so Σ(T − T̄) = 0 and a uniformly warm box
    // does not levitate.
    if b.gravity != 0.0 {
        let a = 1; // gravity is −y
        if t.active[a] {
            for i in 0..n {
                let c = t.coord(i);
                if is_wall(&t, fs, a, c) {
                    continue;
                }
                let (Some(lo), Some(hi)) = (t.cell(step_by(c, a, -1)), t.cell(c)) else {
                    continue;
                };
                let temp = 0.5 * (fs.temperature[lo] + fs.temperature[hi]);
                let alpha = 0.5 * (fs.alpha[lo] + fs.alpha[hi]);
                let force = fs.rho0 * alpha * (temp - fs.t_ref) * b.gravity * vol;
                let q = Momentum::from_si(force * dt).0;
                scratch.delta[a][hi] += q;
                impulse[a] += q;
            }
        }
    }

    // ── Apply. ───────────────────────────────────────────────────────────────
    {
        let f = grid.fluid.as_mut().unwrap();
        let moms: [&mut [i128]; 3] = [f.mom_x, f.mom_y, f.mom_z];
        for (a, m) in moms.into_iter().enumerate() {
            for i in 0..n {
                m[i] += scratch.delta[a][i];
            }
        }
    }

    // ── Project. Make it incompressible. ─────────────────────────────────────
    let pimp = project(grid, b, fs, dt, scratch, cell_mass);
    for a in 0..3 {
        impulse[a] += pimp[a];
    }

    impulse
}


/// Solve `(I − dt·ν·∇²)·u = u₀` for one velocity component, in place.
///
/// Backward Euler on the viscous term, by Gauss–Seidel, exactly as Phase 4b did
/// for heat. It is the *proposal* half of "the solver proposes, the transport
/// disposes": what comes back is a better velocity field, and nothing has moved
/// yet. The caller then runs the ordinary symmetric-flux loop against this
/// field, so conservation stays an identity and does not depend on this
/// converging.
///
/// **Plain Gauss–Seidel, never over-relaxed.** The update is a weighted average
/// of a face's own previous velocity and its neighbours', with positive weights
/// summing to one — a convex combination, which cannot leave the range it
/// started in however few sweeps run. Over-relaxation breaks that, and Phase 4b
/// paid for the lesson with temperatures of 4244 K in a box that began between
/// 400 K and 1200 K.
fn viscous_solve(
    t: &Topo,
    fs: &FluidState,
    a: usize,
    dt: f64,
    dx: f64,
    iters: u32,
    u0: &[f64],
    u: &mut Vec<f64>,
) {
    let n = u0.len();
    u.clear();
    u.extend_from_slice(u0);
    if dt <= 0.0 || dx <= 0.0 {
        return;
    }
    for _ in 0..iters {
        for i in 0..n {
            let c = t.coord(i);
            if is_wall(t, fs, a, c) {
                continue;
            }
            // Assemble the row: `u[i]·(1 + Σβ) − Σβ·u[j] = u0[i]`, with β the
            // face coefficient `dt·ν/dx²`. Solving for `u[i]` gives the convex
            // combination above.
            let mut sum = 0.0;
            let mut diag = 1.0;
            for bx in 0..3 {
                if !t.active[bx] {
                    continue;
                }
                for dir in [-1isize, 1] {
                    let cn = step_by(c, bx, dir);
                    if is_wall(t, fs, a, cn) {
                        // A no-slip wall holds u = 0 there, which stiffens the
                        // row without contributing to the sum — the same
                        // treatment the explicit path gives it.
                        continue;
                    }
                    let Some(j) = t.cell(cn) else { continue };
                    let nu_h = nu_at_face(t, fs, a, c);
                    let nu_t = nu_at_face(t, fs, a, cn);
                    if nu_h <= 0.0 || nu_t <= 0.0 {
                        continue;
                    }
                    let nu = 2.0 * nu_h * nu_t / (nu_h + nu_t).max(1e-300);
                    let beta = dt * nu / (dx * dx);
                    diag += beta;
                    sum += beta * u[j];
                }
            }
            if diag > 0.0 {
                u[i] = (u0[i] + sum) / diag;
            }
        }
    }
}

#[inline]
fn nu_at_face(t: &Topo, fs: &FluidState, a: usize, c: [isize; 3]) -> f64 {
    match (t.cell(step_by(c, a, -1)), t.cell(c)) {
        (Some(l), Some(h)) => 0.5 * (fs.nu[l] + fs.nu[h]),
        _ => 0.0,
    }
}

/// Reusable working memory for the fluid step.
#[derive(Default)]
pub struct FluidScratch {
    delta: [Vec<i128>; 3],
    /// Face velocities before the viscous solve, and the field it proposes.
    u0: [Vec<f64>; 3],
    u_new: [Vec<f64>; 3],
    div: Vec<f64>,
    /// The Poisson stencil, precomputed once per substep.
    ///
    /// The pressure solve sweeps the grid sixty times; recomputing "which of my
    /// neighbours are real" inside that loop meant three integer divisions and
    /// six wall tests per cell per sweep, and it was three-quarters of the total
    /// runtime of the entire engine. Working out the answer once and looking it
    /// up sixty times is the whole optimisation, and it is a factor of ten.
    stencil: Vec<[i32; 6]>,
    stencil_n: Vec<f64>,
    coords: Vec<[isize; 3]>,
    residual: f64,
    max_div: f64,
}

impl FluidScratch {
    pub fn new() -> Self {
        Self::default()
    }
    fn resize(&mut self, n: usize) {
        for d in self.delta.iter_mut() {
            d.clear();
            d.resize(n, 0);
        }
        self.div.clear();
        self.div.resize(n, 0.0);
        self.stencil.clear();
        self.stencil.resize(n, [-1; 6]);
        self.stencil_n.clear();
        self.stencil_n.resize(n, 0.0);
        self.coords.clear();
        self.coords.resize(n, [0; 3]);
    }
    /// How badly the last pressure solve failed to converge. Reported, never
    /// chased — see the note on `projection_iters`.
    pub fn residual(&self) -> f64 {
        self.residual
    }
    /// The worst |∇·u| left in the field. This is the honest measure of whether
    /// the fluid is actually incompressible or merely nearly so.
    pub fn max_divergence(&self) -> f64 {
        self.max_div
    }
}

/// The pressure projection (Chorin, 1968).
///
/// # The dominant cost of the whole engine, and why it cannot simply be cut
///
/// Measured on a 24x12 protocell convection cell: twenty ticks cost 4.09 s at
/// sixty sweeps and 0.20 s at fifteen. Twentyfold, for a fourfold cut — the
/// projection is not *a* cost, it is very nearly *the* cost.
///
/// And the cheap win is not available. Measured divergence, as a fraction of the
/// velocity scale:
///
/// ```text
///     15 sweeps    200.5%      grossly compressible
///     30 sweeps     91.0%
///     60 sweeps      1.4%      honest
/// ```
///
/// Convergence is sharply nonlinear and sixty sweeps is near the edge of it, not
/// comfortable padding. An under-converged projection leaves matter accumulating
/// where nothing put it, and the buoyancy that drifted density produces looks
/// *exactly* like convection — the artefact that produced Phase 2b's fictitious
/// stirring at half the critical Rayleigh number, with the perturbation set to
/// zero.
///
/// So making fine grids affordable needs a better **method**, not fewer sweeps.
/// Single-grid Gauss-Seidel on a Laplace problem converges at a rate set by the
/// *grid* — roughly `cos^2(pi/n)` per sweep — which is the same finding Phase 5b
/// wrote down about deep-time diffusion, where the fix was `2*n^2` sweeps and the
/// note that multigrid is the real answer. It is the real answer here too, and it
/// would serve conduction, species transport and this projection together.
///
/// After advection, viscosity and buoyancy, the velocity field is *wrong*: it has
/// divergence, which means fluid is appearing and disappearing. Find the pressure
/// field whose gradient removes exactly that divergence, and subtract it.
///
/// ```text
/// ∇²p = (ρ₀/dt) ∇·u*        then      u = u* − (dt/ρ₀) ∇p
/// ```
///
/// Solved with red-black SOR. Red-black because the two colours never touch, so
/// every cell of a colour can be updated against a fixed set of neighbours — which
/// is both deterministic and, when it matters, parallel. SOR because plain
/// Gauss–Seidel would need an order of magnitude more sweeps.
///
/// With walls and periodicity everywhere and no open boundary, the system is
/// **singular**: p is defined only up to a constant, and the right-hand side must
/// sum to zero or there is no solution at all. Both are handled below, and
/// forgetting either produces a solver that diverges slowly and mysteriously.
fn project(
    grid: &mut GridView,
    b: &Boundary,
    fs: &FluidState,
    dt: f64,
    scratch: &mut FluidScratch,
    cell_mass: f64,
) -> [i128; 3] {
    let t = Topo::new(grid.spec, b);
    let n = grid.spec.cells();
    let dx = grid.spec.cell_m;
    let vol = grid.spec.cell_volume();

    // The stencil, once. A face that is a wall imposes ∂p/∂n = 0 — its ghost value
    // equals our own, so it contributes nothing and does not count. A solid
    // neighbour is a wall. That single sentence is the entire immersed-boundary
    // treatment, and it is why a convecting magma layer under a frozen lid needs
    // no special-case code at all.
    for i in 0..n {
        scratch.coords[i] = t.coord(i);
    }
    for i in 0..n {
        let c = scratch.coords[i];
        let mut k = 0usize;
        let mut st = [-1i32; 6];
        if !fs.solid[i] {
            for a in 0..3 {
                if !t.active[a] {
                    continue;
                }
                for d in [-1isize, 1] {
                    let face = if d > 0 { step_by(c, a, 1) } else { c };
                    if is_wall(&t, fs, a, face) {
                        continue;
                    }
                    if let Some(j) = t.cell(step_by(c, a, d)) {
                        st[k] = j as i32;
                        k += 1;
                    }
                }
            }
        }
        scratch.stencil[i] = st;
        scratch.stencil_n[i] = k as f64;
    }

    // ∇·u*
    let (mx, my, mz) = {
        let f = grid.fluid.as_ref().unwrap();
        (f.mom_x.to_vec(), f.mom_y.to_vec(), f.mom_z.to_vec())
    };
    let mom_ro: [&[i128]; 3] = [&mx, &my, &mz];

    let mut rhs_sum = 0.0;
    let mut fluid_cells = 0.0;
    for i in 0..n {
        if fs.solid[i] {
            scratch.div[i] = 0.0;
            continue;
        }
        let c = t.coord(i);
        let mut d = 0.0;
        for a in 0..3 {
            if !t.active[a] {
                continue;
            }
            d += (u_face(&t, fs, &mom_ro, a, step_by(c, a, 1))
                - u_face(&t, fs, &mom_ro, a, c))
                / dx;
        }
        scratch.div[i] = d;
        rhs_sum += d;
        fluid_cells += 1.0;
    }

    // Compatibility. An all-Neumann Poisson problem has a solution only if the
    // source integrates to zero — you cannot inflate a sealed box. Any residual
    // here is accumulated numerical error, and subtracting its mean is the
    // standard, and necessary, cure.
    let mean = if fluid_cells > 0.0 { rhs_sum / fluid_cells } else { 0.0 };
    let scale = fs.rho0 / dt;
    for i in 0..n {
        if !fs.solid[i] {
            scratch.div[i] = (scratch.div[i] - mean) * scale;
        }
    }

    // Red-black SOR, warm-started from the previous step's pressure.
    let p: &mut [f64] = grid.fluid.as_mut().unwrap().pressure;
    let dx2 = dx * dx;
    // Zero means "derive it" — see `Boundary::sor_omega`. The relevant length is
    // the *smallest* active dimension, because convergence is set by the slowest
    // direction to communicate across.
    let omega = if b.sor_omega > 0.0 {
        b.sor_omega
    } else {
        let n = [grid.spec.nx, grid.spec.ny, grid.spec.nz]
            .into_iter()
            .filter(|&d| d > 1)
            .min()
            .unwrap_or(1);
        crate::grid::derived_sor_omega(n)
    };

    for _ in 0..b.projection_iters {
        // Red-black: the two colours never touch, so a whole colour can be
        // updated against a frozen set of neighbours. Deterministic, and — when
        // it matters — parallel for free.
        for colour in 0..2 {
            for i in 0..n {
                let cnt = scratch.stencil_n[i];
                if cnt == 0.0 {
                    continue;
                }
                let c = scratch.coords[i];
                if ((c[0] + c[1] + c[2]) & 1) as usize != colour {
                    continue;
                }
                let st = &scratch.stencil[i];
                let mut sum = 0.0;
                for k in 0..(cnt as usize) {
                    sum += p[st[k] as usize];
                }
                let target = (sum - dx2 * scratch.div[i]) / cnt;
                p[i] = (1.0 - omega) * p[i] + omega * target;
            }
        }
    }

    // Pin the gauge. p is determined only up to a constant; letting it wander
    // does no physical harm but it drifts without bound, and an unbounded number
    // in the world hash is a slow-motion overflow.
    let psum: f64 = (0..n).filter(|&i| !fs.solid[i]).map(|i| p[i]).sum();
    if fluid_cells > 0.0 {
        let pm = psum / fluid_cells;
        for i in 0..n {
            if !fs.solid[i] {
                p[i] -= pm;
            }
        }
    }

    // Residual — how well did it actually solve? Reported, not chased.
    let mut res: f64 = 0.0;
    for i in 0..n {
        let cnt = scratch.stencil_n[i] as usize;
        if cnt == 0 {
            continue;
        }
        let st = &scratch.stencil[i];
        let mut lap = 0.0;
        for k in 0..cnt {
            lap += (p[st[k] as usize] - p[i]) / dx2;
        }
        res = res.max((lap - scratch.div[i]).abs());
    }
    scratch.residual = res;

    // Correct the momentum: subtract the pressure gradient.
    let mut impulse = [0i128; 3];
    let pv: Vec<f64> = p.to_vec();
    let _ = &p; // borrow released at end of scope; pv holds the copy we need
    let f = grid.fluid.as_mut().unwrap();
    let moms: [&mut [i128]; 3] = [f.mom_x, f.mom_y, f.mom_z];
    for (a, m) in moms.into_iter().enumerate() {
        if !t.active[a] {
            continue;
        }
        for i in 0..n {
            let c = t.coord(i);
            if is_wall(&t, fs, a, c) {
                continue;
            }
            let (Some(lo), Some(hi)) = (t.cell(step_by(c, a, -1)), t.cell(c)) else {
                continue;
            };
            let grad = (pv[hi] - pv[lo]) / dx;
            let q = Momentum::from_si(-grad * vol * dt).0;
            m[i] += q;
            // Summed over the domain this is the net force the boundary exerts on
            // the fluid. Not an internal exchange, so it goes in the books.
            impulse[a] += q;
        }
    }

    // How incompressible did we actually manage to be?
    let (mx, my, mz) = {
        let f = grid.fluid.as_ref().unwrap();
        (f.mom_x.to_vec(), f.mom_y.to_vec(), f.mom_z.to_vec())
    };
    let mom_ro: [&[i128]; 3] = [&mx, &my, &mz];
    let mut md: f64 = 0.0;
    for i in 0..n {
        if fs.solid[i] {
            continue;
        }
        let c = t.coord(i);
        let mut d = 0.0;
        for a in 0..3 {
            if !t.active[a] {
                continue;
            }
            d += (u_face(&t, fs, &mom_ro, a, step_by(c, a, 1))
                - u_face(&t, fs, &mom_ro, a, c))
                / dx;
        }
        md = md.max(d.abs());
    }
    scratch.max_div = md;
    let _ = cell_mass;

    impulse
}

/// The velocity field, for advecting scalars and for anyone who wants to look.
pub fn velocities(grid: &GridView, b: &Boundary, fs: &FluidState) -> Vec<[f64; 3]> {
    let t = Topo::new(grid.spec, b);
    let n = grid.spec.cells();
    let Some(f) = grid.fluid.as_ref() else { return vec![[0.0; 3]; n] };
    let mom_ro: [&[i128]; 3] = [f.mom_x, f.mom_y, f.mom_z];
    let mut out = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let c = t.coord(i);
        for a in 0..3 {
            // Cell-centred velocity: the mean of the two faces. For looking at,
            // not for solving with.
            out[i][a] = 0.5
                * (u_face(&t, fs, &mom_ro, a, c) + u_face(&t, fs, &mom_ro, a, step_by(c, a, 1)));
        }
    }
    out
}

/// Face velocities, indexed as the −a face of each cell. What the scalar
/// advection actually uses — the divergence-free field, not the interpolated one.
pub fn face_velocities(grid: &GridView, b: &Boundary, fs: &FluidState) -> [Vec<f64>; 3] {
    let t = Topo::new(grid.spec, b);
    let n = grid.spec.cells();
    let mut out = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    let Some(f) = grid.fluid.as_ref() else { return out };
    let mom_ro: [&[i128]; 3] = [f.mom_x, f.mom_y, f.mom_z];
    for a in 0..3 {
        for i in 0..n {
            out[a][i] = u_face(&t, fs, &mom_ro, a, t.coord(i));
        }
    }
    out
}

/// Largest |u| anywhere. Sets the advective timestep limit: a parcel may not
/// cross more than one cell per step, or the scheme is solving a different
/// problem from the one it was handed.
pub fn max_speed(grid: &GridView, b: &Boundary, fs: &FluidState) -> f64 {
    face_velocities(grid, b, fs)
        .iter()
        .flat_map(|v| v.iter())
        .fold(0.0f64, |m, &x| m.max(x.abs()))
}

/// Total kinetic energy, J. Not a conserved quantity — viscosity eats it — but
/// the thing you watch to see whether a convective instability is growing or
/// dying.
pub fn kinetic_energy(grid: &GridView, b: &Boundary, fs: &FluidState) -> f64 {
    if fs.rho0 <= 0.0 {
        return 0.0;
    }
    let m = fs.rho0 * grid.spec.cell_volume();
    face_velocities(grid, b, fs)
        .iter()
        .flat_map(|v| v.iter())
        .map(|&u| 0.5 * m * u * u)
        .sum()
}

/// Which cells are rigid. For the renderer, and for anyone who wants to see where
/// the magma stops and the crust begins.
pub fn solid_mask(fs: &FluidState) -> &[bool] {
    &fs.solid
}

#[allow(dead_code)]
fn _phase_is_used(_: Phase) {}
