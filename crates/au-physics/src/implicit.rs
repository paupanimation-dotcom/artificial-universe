//! Implicit conduction — the phase that unlocks time.
//!
//! # The crack this closes
//!
//! Every earlier phase has been honest about a limit it could not pass: an
//! explicit diffusion scheme is stable only while
//!
//! ```text
//!     dt  ≤  dx² / (2·ndim·α)
//! ```
//!
//! and that `dx²` is fatal. Halve the cell size for detail and the timestep must
//! quarter; the cost of a fixed span of simulated time rises as the *fourth* power
//! of resolution. Phase 1's deep-time accelerator turned out to be largely inert
//! against real physics for exactly this reason. Weather needs hours, geology
//! needs epochs, and evolution needs generations — none of them reachable while
//! the timestep is chained to the square of the grid spacing.
//!
//! Backward Euler breaks the chain. Instead of computing the new temperature from
//! the old one:
//!
//! ```text
//!     explicit:   T^{n+1} = T^n + dt·α·∇²T^n         — stable only if dt is small
//!     implicit:   T^{n+1} = T^n + dt·α·∇²T^{n+1}     — unconditionally stable
//! ```
//!
//! the implicit form defines the new field in terms of itself, which turns one
//! explicit sweep into a sparse linear system, `(I − dt·α·∇²)·T^{n+1} = T^n`. That
//! system is diagonally dominant and nearest-neighbour — the same shape as the
//! pressure projection Phase 2b already solves — so the red-black SOR built there
//! solves this too. The reward is that no timestep is unstable. A step ten
//! thousand times past the explicit limit still produces a bounded, monotone,
//! physically sane field; it is merely *less accurate*, which is a trade a caller
//! can make knowingly.
//!
//! # The hard part: an approximate solver inside an exact engine
//!
//! Everything in this project rests on energy being an exact `i128`, conserved to
//! the microjoule. An implicit solver is iterative floating point that converges
//! to a tolerance and never lands exactly. Bolted in naively, it would put a
//! rounding leak at the centre of the foundation — and a leak is precisely what a
//! later evolutionary optimiser learns to farm.
//!
//! The resolution is to give the solver and the transport **different jobs**:
//!
//!   * The **solver** decides *how much energy should move*. Floating point,
//!     iterative, approximate, allowed to be imperfect.
//!   * The **transport** moves it: for each face, one integer, subtracted from one
//!     side and added to the other — `delta[a] -= de; delta[b] += de`.
//!
//! So **conservation does not depend on the solve converging.** An
//! under-converged solve puts heat in slightly the wrong place — wrong physics,
//! honestly wrong, and visible in the reported residual. It cannot lose a joule,
//! because losing one would require the two halves of a symmetric integer transfer
//! to disagree, and they are the same integer.
//!
//! This is the third time the project has needed the same idea. Phase 2 made
//! conduction conservative by construction with symmetric flux; Phase 3 fixed the
//! combustion dynamic-range bug by quantising the *sum* and never the summand;
//! here, approximation is again allowed to make the answer wrong but never allowed
//! to make it non-conservative.
//!
//! The corollary is a trap worth naming, because it is the obvious implementation
//! and it is wrong: do **not** apply the energy change as `C·(T^{n+1} − T^n)` per
//! cell. That is conservative only if the solve converged perfectly, so it silently
//! couples the foundation's central invariant to a tolerance parameter. Apply face
//! fluxes instead — conservative always, converged or not.
//!
//! # Scope, and what stays explicit
//!
//! Interior conduction goes implicit here. That is the term with the `dx²` limit,
//! and therefore the one worth the machinery. Deliberately left explicit:
//!
//!   * **Radiation** — nonlinear in `T⁴`; making it implicit needs Newton
//!     iteration. Its stability bound scales as `T³`, not `1/dx²`, so it does not
//!     worsen when the grid is refined, and the existing limiter already covers it.
//!   * **Advection** — bounded by `dt < dx/u`, linear in `dx`. Uncomfortable, not
//!     a cliff.
//!   * **Viscosity** — has the same `dx²` problem as heat and deserves the same
//!     cure; it is the next candidate, and the machinery here is written to be
//!     reusable for it.

use crate::grid::GridSpec;
use crate::quantity::Energy;
use crate::state::CellState;

/// Working memory for the implicit solve. Reused across steps so a long run does
/// no allocation; holds no state that survives a call, only scratch.
#[derive(Default)]
pub struct ImplicitScratch {
    /// The iterate — the temperature field being solved for.
    t: Vec<f64>,
    /// `T^n`, the field at the start of the step.
    t0: Vec<f64>,
    /// Total heat capacity per cell, J/K.
    cap: Vec<f64>,
    /// Neighbour cell indices, up to six per cell.
    nb: Vec<[u32; 6]>,
    /// Conductance to each neighbour, W/K (harmonic mean, as in explicit).
    g: Vec<[f64; 6]>,
    /// How many neighbours this cell actually has.
    nbn: Vec<u8>,
    /// Cells that take part: finite heat capacity and some conductivity.
    active: Vec<bool>,
    /// Parity of each cell, for red-black ordering.
    parity: Vec<u8>,

    /// Largest absolute residual of the linear system after the final sweep, K.
    /// Reported, never chased — see the note on fixed iteration counts below.
    pub residual: f64,
}

impl ImplicitScratch {
    pub fn new() -> Self {
        Self::default()
    }

    fn resize(&mut self, n: usize) {
        self.t.clear();
        self.t.resize(n, 0.0);
        self.t0.clear();
        self.t0.resize(n, 0.0);
        self.cap.clear();
        self.cap.resize(n, 0.0);
        self.nb.clear();
        self.nb.resize(n, [0; 6]);
        self.g.clear();
        self.g.resize(n, [0.0; 6]);
        self.nbn.clear();
        self.nbn.resize(n, 0);
        self.active.clear();
        self.active.resize(n, false);
        self.parity.clear();
        self.parity.resize(n, 0);
    }
}

/// Advance conduction implicitly over `dt`, writing the resulting energy changes
/// into `delta` as exact integer symmetric transfers.
///
/// `iters` is a **fixed** sweep count, not a convergence criterion. A loop that
/// stops "when converged" makes the amount of work — and therefore the answer —
/// depend on the state, which would make the simulation's output a function of its
/// own numerical history in a way nothing else in this engine is. The residual is
/// measured and reported instead, so a caller who is not converging finds out
/// rather than silently paying for it. (The pressure projection made the same
/// choice for the same reason.)
///
/// `omega` is the SOR over-relaxation factor; 1.0 is plain Gauss–Seidel.
pub fn conduct_implicit(
    spec: GridSpec,
    states: &[CellState],
    area: f64,
    dx: f64,
    dt: f64,
    iters: u32,
    omega: f64,
    scratch: &mut ImplicitScratch,
    delta: &mut [i128],
) {
    let n = spec.cells();
    if n == 0 || dt <= 0.0 {
        return;
    }
    scratch.resize(n);

    // ── Assemble. Conductance across a face is the harmonic mean of the two
    //    conductivities — two slabs in series add their *resistances* — which is
    //    the identical rule the explicit path uses. That shared rule is why the
    //    two schemes must agree in the limit of small dt, and a test holds them
    //    to it.
    for z in 0..spec.nz {
        for y in 0..spec.ny {
            for x in 0..spec.nx {
                let i = spec.idx(x, y, z);
                let s = states[i];
                scratch.t0[i] = s.temperature;
                scratch.t[i] = s.temperature; // warm start from the current field
                scratch.cap[i] = s.heat_capacity;
                scratch.parity[i] = ((x + y + z) & 1) as u8;
                scratch.active[i] = s.heat_capacity > 0.0 && s.conductivity > 0.0;
            }
        }
    }

    for z in 0..spec.nz {
        for y in 0..spec.ny {
            for x in 0..spec.nx {
                let i = spec.idx(x, y, z);
                if !scratch.active[i] {
                    continue;
                }
                let ka = states[i].conductivity;
                let mut count = 0usize;
                let push = |j: usize, count: &mut usize, nb: &mut [u32; 6], g: &mut [f64; 6]| {
                    if !scratch.active[j] {
                        return;
                    }
                    let kb = states[j].conductivity;
                    let k = 2.0 * ka * kb / (ka + kb);
                    nb[*count] = j as u32;
                    g[*count] = k * area / dx;
                    *count += 1;
                };
                let mut nb = [0u32; 6];
                let mut g = [0.0f64; 6];
                if x > 0 {
                    push(spec.idx(x - 1, y, z), &mut count, &mut nb, &mut g);
                }
                if x + 1 < spec.nx {
                    push(spec.idx(x + 1, y, z), &mut count, &mut nb, &mut g);
                }
                if y > 0 {
                    push(spec.idx(x, y - 1, z), &mut count, &mut nb, &mut g);
                }
                if y + 1 < spec.ny {
                    push(spec.idx(x, y + 1, z), &mut count, &mut nb, &mut g);
                }
                if z > 0 {
                    push(spec.idx(x, y, z - 1), &mut count, &mut nb, &mut g);
                }
                if z + 1 < spec.nz {
                    push(spec.idx(x, y, z + 1), &mut count, &mut nb, &mut g);
                }
                scratch.nb[i] = nb;
                scratch.g[i] = g;
                scratch.nbn[i] = count as u8;
            }
        }
    }

    // ── Solve. Per cell the backward-Euler row is
    //
    //      C·(T_i − T⁰_i)/dt  =  Σ_j G_ij·(T_j − T_i)
    //
    //    which rearranges to the update below. The denominator exceeds the sum of
    //    the off-diagonal weights for every positive dt, so the matrix is strictly
    //    diagonally dominant and this iteration converges for *any* timestep —
    //    that unconditional property is the entire reason for the phase.
    for _ in 0..iters {
        for colour in 0..2u8 {
            for i in 0..n {
                if !scratch.active[i] || scratch.parity[i] != colour {
                    continue;
                }
                let cnt = scratch.nbn[i] as usize;
                if cnt == 0 {
                    continue;
                }
                let cap = scratch.cap[i];
                if cap <= 0.0 {
                    continue;
                }
                let mut gsum = 0.0;
                let mut gt = 0.0;
                for k in 0..cnt {
                    let gk = scratch.g[i][k];
                    gsum += gk;
                    gt += gk * scratch.t[scratch.nb[i][k] as usize];
                }
                let f = dt / cap;
                let target = (scratch.t0[i] + f * gt) / (1.0 + f * gsum);
                scratch.t[i] = (1.0 - omega) * scratch.t[i] + omega * target;
            }
        }
    }

    // ── Report how well it actually solved, in kelvin. Not acted on: a caller
    //    that sees a large residual should raise the iteration count or shorten
    //    the step, and the decision is theirs to make explicitly.
    let mut res: f64 = 0.0;
    for i in 0..n {
        if !scratch.active[i] {
            continue;
        }
        let cnt = scratch.nbn[i] as usize;
        let cap = scratch.cap[i];
        if cnt == 0 || cap <= 0.0 {
            continue;
        }
        let mut acc = 0.0;
        for k in 0..cnt {
            acc += scratch.g[i][k] * (scratch.t[scratch.nb[i][k] as usize] - scratch.t[i]);
        }
        let r = (scratch.t[i] - scratch.t0[i]) - (dt / cap) * acc;
        if r.abs() > res {
            res = r.abs();
        }
    }
    scratch.residual = res;

    // ── Transport. **This** is where conservation lives, and it lives here
    //    whether or not the solve above succeeded. Each face is visited once, in
    //    the +x/+y/+z direction only, so every pair is handled exactly one time;
    //    the energy that crosses is a single integer removed from one cell and
    //    added to the other. No tolerance, no accumulation, no drift.
    let face = |a: usize, bi: usize, delta: &mut [i128]| {
        if !scratch.active[a] || !scratch.active[bi] {
            return;
        }
        let (ka, kb) = (states[a].conductivity, states[bi].conductivity);
        let k = 2.0 * ka * kb / (ka + kb);
        let dtemp = scratch.t[a] - scratch.t[bi];
        if dtemp == 0.0 {
            return;
        }
        let de = Energy::from_joules(k * area / dx * dtemp * dt).0;
        delta[a] -= de;
        delta[bi] += de;
    };

    for z in 0..spec.nz {
        for y in 0..spec.ny {
            for x in 0..spec.nx {
                let i = spec.idx(x, y, z);
                if x + 1 < spec.nx {
                    face(i, spec.idx(x + 1, y, z), delta);
                }
                if y + 1 < spec.ny {
                    face(i, spec.idx(x, y + 1, z), delta);
                }
                if z + 1 < spec.nz {
                    face(i, spec.idx(x, y, z + 1), delta);
                }
            }
        }
    }
}
