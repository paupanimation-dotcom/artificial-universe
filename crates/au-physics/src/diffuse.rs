//! Diffusion of a conserved integer count — the physics under Phase 5b.
//!
//! # What this is for
//!
//! Phase 5a gave the world molecules it was never told about; this module gives
//! them legs. A discovered species is a population count per cell, and Fick's
//! law says counts flow down their own gradient: through a face between two
//! cells of equal volume, the transfer rate is
//!
//! ```text
//!     dN/dt = (D / dx²) · (N_a − N_b)
//! ```
//!
//! — exactly conduction's shape with concentration in temperature's chair, which
//! is why this file is `implicit.rs`'s deliberate twin. It shares the two
//! commitments that made conduction trustworthy:
//!
//!   * **The solver proposes, the transport disposes.** Whatever field the
//!     (possibly approximate, possibly implicit) solve produces, matter moves
//!     only as one integer per face, subtracted from one side and added to the
//!     other. Conservation never depends on convergence. Third verse, same as
//!     the first: symmetric flux (Phase 2), quantise the sum (Phase 3), and now
//!     this.
//!
//!   * **Unconditional stability is available, because deep time demands it.**
//!     Explicit diffusion carries the same `dt ≤ dx²/(2·d·D)` cliff conduction
//!     had, and Phase 4b's cure — backward Euler, plain Gauss–Seidel, fixed
//!     sweeps, ω = 1 so the update is a convex combination and the maximum
//!     principle survives any under-convergence — transplants unchanged. (This
//!     is the "written to be reusable" promise from Phase 4b being kept; a
//!     shared generic core between the two is a refactor that waits for a third
//!     user, per the rule of three.)
//!
//! # Counts, not concentrations
//!
//! Cells in one grid share a volume, so `N_a − N_b` and `c_a − c_b` differ by a
//! constant factor that folds into the rate. Working in raw counts keeps the
//! transported quantity the same exact integer the chemistry columns store —
//! no unit conversion, no rounding seam between storage and transport.
//!
//! # The mobility mask, and what it buys
//!
//! Species here are **solutes**: they live dissolved in a cell's host material
//! and move only where that host is fluid. A face conducts species only if
//! *both* cells are mobile (liquid or gas host, as derived by the same
//! `au_physics::derive` everything else trusts). Solids and vacuum are walls.
//!
//! That one rule is not a technicality — it is an emergence lever. A cooling
//! ocean that crosses its freezing point *traps its chemistry in place*: the
//! spatial pattern of dissolved species at the moment of freezing becomes a
//! frozen stratigraphy, and a later thaw releases it. Nobody will script that;
//! it falls out of "solutes need a fluid host" meeting "phase is derived from
//! energy".
//!
//! # The quantisation deadband, stated honestly
//!
//! Transfers round to the nearest integer, so a face whose ideal flux is below
//! half a molecule per step moves nothing: gradients shallower than about
//! `1/(2k)` counts per face persist indefinitely. At the population scales the
//! engine runs (10⁶ and up) this is invisible; at tens of molecules it is a
//! real floor. The honest microphysics at that scale is stochastic — individual
//! Brownian hops — and doing it honestly means drawing from the derived RNG,
//! which would make diffusion the first *random* physics in the engine. That is
//! a real decision with real consequences for the determinism story, so it is
//! deferred and documented rather than smuggled in via rounding tricks.
//!
//! # How many sweeps deep time actually needs
//!
//! This was learned the hard way, so it is written down. At large `k` the
//! implicit system approaches a pure Laplace problem, and single-grid
//! Gauss–Seidel's convergence factor there is set by the *grid*, not by `k`:
//! roughly `cos²(π/n)` per sweep along the longest dimension `n`. Information
//! crawls one cell per sweep, so driving the residual down takes **O(n²)
//! sweeps** — 24 sweeps that comfortably converge a small-k step leave a
//! k ≈ 10³ system on a 41-cell tube with rows violated by ~10⁶ counts.
//!
//! And an under-converged solve here is worse than a slow one, because the
//! transport multiplies *face differences* by `k`: solver error ε becomes
//! transported error `k·ε`. In the first deep-time runs that noise manufactured
//! a red-black checkerboard in the populations, which the warm start then
//! re-imprinted on the next stride's solve — a self-perpetuating artifact,
//! ends starved to exactly zero, diagnosed with a seven-cell probe. The fix is
//! not cleverness but honesty about the method: sweeps must scale as `2·n²`
//! (`sweeps_for_deep_time` below), after which the same 41-tube mixes to the
//! *exact* uniform count with residual 0.0. Multigrid would make the cost
//! O(n) instead of O(n³) total and is noted as future work; at chunk-scale
//! grids, paying n² is cheap and simple.
//!
//! Conduction's solver (`implicit.rs`) shares this trait — its validated
//! regimes sit below the threshold where it bites, and its transport divides
//! by heat capacity, which blunts the amplification. Unifying both under one
//! multigrid core is the natural future; until then this note is the fence.
//!
//! # The donor's-pocket limiter
//!
//! An under-converged implicit solve can propose face fluxes that overdraw a
//! cell — neighbouring solved values still differ by O(curvature), and
//! multiplied by an enormous `k` that difference can exceed everything the
//! donor holds. The sum would still balance (symmetric transfers cannot lose
//! matter), but a *negative population* is impossible physics, and the first
//! test that ran one sweep at k = 300 produced exactly that.
//!
//! The guard is the finite-volume classic: **no face may export more than the
//! donor's holdings divided by the donor's open-face count.** Then even if
//! every face draws its maximum, the cell empties to zero and no further. The
//! cap reads only the frozen pre-step field and a fixed degree, so it is
//! order-independent and cannot break symmetry or conservation — and whenever
//! the solve is honest (converged fluxes are small), it never engages. A
//! guard, not a scheme: it makes wrong physics safely wrong instead of
//! impossibly wrong.

use crate::grid::GridSpec;

/// Working memory for the implicit solve. Reused across calls; no state
/// survives a call.
#[derive(Default)]
pub struct DiffuseScratch {
    /// The iterate — the count field being solved for, as f64.
    n: Vec<f64>,
    /// The field at the start of the step.
    n0: Vec<f64>,
    /// Red-black parity per cell.
    parity: Vec<u8>,
    /// The pre-step field, frozen: what each donor actually holds.
    frozen: Vec<i128>,
    /// Open (mobile↔mobile) face count per cell — the divisor in the limiter.
    degree: Vec<u8>,
    /// Largest |residual| of the linear system after the final sweep, in
    /// counts. Reported, never chased — fixed sweep counts keep the work (and
    /// therefore the answer) a pure function of the state.
    pub residual: f64,
}

impl DiffuseScratch {
    pub fn new() -> Self {
        Self::default()
    }
    fn resize(&mut self, n: usize) {
        self.n.clear();
        self.n.resize(n, 0.0);
        self.n0.clear();
        self.n0.resize(n, 0.0);
        self.parity.clear();
        self.parity.resize(n, 0);
        self.frozen.clear();
        self.frozen.resize(n, 0);
        self.degree.clear();
        self.degree.resize(n, 0);
    }
}

/// Visit every interior face once (+x, +y, +z), calling `f(a, b)` on the two
/// cell indices. One visit per face is what makes symmetric transfer exact.
#[inline]
fn for_each_face(spec: GridSpec, mut f: impl FnMut(usize, usize)) {
    for z in 0..spec.nz {
        for y in 0..spec.ny {
            for x in 0..spec.nx {
                let i = spec.idx(x, y, z);
                if x + 1 < spec.nx {
                    f(i, spec.idx(x + 1, y, z));
                }
                if y + 1 < spec.ny {
                    f(i, spec.idx(x, y + 1, z));
                }
                if z + 1 < spec.nz {
                    f(i, spec.idx(x, y, z + 1));
                }
            }
        }
    }
}

/// The sweep count deep time needs: `max(floor, 2·L²)` where `L` is the
/// longest grid dimension. See the module docs ("How many sweeps…") for the
/// convergence argument; the short version is that Gauss–Seidel moves
/// information one cell per sweep and its Laplace-limit convergence factor is
/// `cos²(π/L)`, so the sweep budget must grow with the square of the span.
pub fn sweeps_for_deep_time(spec: GridSpec, floor: u32) -> u32 {
    let l = spec.nx.max(spec.ny).max(spec.nz) as u32;
    floor.max(2 * l * l)
}

/// The permeability of a cell: 1.0 when nothing is in the way.
///
/// A membrane does not close a face, it *slows* it — so permeability multiplies
/// the rate rather than clearing the mobility mask. Crossing from one cell to
/// another passes through both barriers, so a face carries the product of the
/// two: `k_ij = k · p_i · p_j`. Symmetric by construction, which is what keeps
/// the transfer exactly antisymmetric and conservation untouched.
#[inline]
fn perm_of(perm: Option<&[f64]>, i: usize) -> f64 {
    match perm {
        Some(p) => p[i],
        None => 1.0,
    }
}

/// Fill `frozen` and `degree` for the current field and mask.
fn prepare(spec: GridSpec, field: &[i128], mobile: &[bool], scratch: &mut DiffuseScratch) {
    scratch.resize(spec.cells());
    scratch.frozen.copy_from_slice(field);
    for_each_face(spec, |a, b| {
        if mobile[a] && mobile[b] {
            scratch.degree[a] += 1;
            scratch.degree[b] += 1;
        }
    });
}

/// The limited symmetric transfer for one face: the raw proposal, capped by the
/// donor's pocket. Reads only frozen state — order-independent by construction.
#[inline]
fn transfer(raw: f64, a: usize, b: usize, scratch: &DiffuseScratch, field: &mut [i128]) {
    let mut de = raw.round() as i128;
    if de > 0 {
        let cap = scratch.frozen[a] / scratch.degree[a].max(1) as i128;
        de = de.min(cap.max(0));
    } else if de < 0 {
        let cap = scratch.frozen[b] / scratch.degree[b].max(1) as i128;
        de = de.max(-cap.max(0));
    }
    field[a] -= de;
    field[b] += de;
}

/// One explicit diffusion step at rate `k = D·dt/dx²` (dimensionless).
///
/// Stability requires `k ≤ 1/(2·ndim)`; the *caller* substeps to honour that,
/// because the caller knows its dt policy. Flux across each open face is
/// `round(k · (N_a − N_b))`, applied as one integer both ways — so the sum over
/// the grid is invariant to the last count, whatever `k` was.
pub fn diffuse_explicit(
    spec: GridSpec,
    field: &mut [i128],
    mobile: &[bool],
    k: f64,
    scratch: &mut DiffuseScratch,
) {
    diffuse_explicit_perm(spec, field, mobile, None, k, scratch)
}

/// One explicit diffusion step with a per-cell permeability barrier.
///
/// `perm[i]` scales how freely cell *i* exchanges with its neighbours; 1.0 is an
/// open cell and 0.0 a sealed one. This is what a membrane does: it does not
/// remove a cell from the world, it makes the world reach it more slowly.
pub fn diffuse_explicit_perm(
    spec: GridSpec,
    field: &mut [i128],
    mobile: &[bool],
    perm: Option<&[f64]>,
    k: f64,
    scratch: &mut DiffuseScratch,
) {
    debug_assert_eq!(field.len(), spec.cells());
    // Fluxes are computed against the frozen pre-step field, then applied —
    // otherwise a sweep order would leak into the answer.
    prepare(spec, field, mobile, scratch);
    for_each_face(spec, |a, b| {
        if !mobile[a] || !mobile[b] {
            return;
        }
        let kf = k * perm_of(perm, a) * perm_of(perm, b);
        let raw = kf * (scratch.frozen[a] - scratch.frozen[b]) as f64;
        transfer(raw, a, b, scratch, field);
    });
}

/// Advance diffusion implicitly (backward Euler) over a step whose rate is
/// `k = D·dt/dx²`, with **no** stability ceiling on `k`.
///
/// The row for cell *i* is `N_i − N⁰_i = k·Σ_j (N_j − N_i)` over open faces —
/// strictly diagonally dominant for every positive `k`, so Gauss–Seidel
/// converges for *any* step size. `omega` should stay at 1.0: the update is
/// then a convex combination of the cell's own start value and its neighbours,
/// which keeps every iterate inside the initial range no matter how few sweeps
/// run — the same maximum-principle argument, and the same 4244 K lesson,
/// as conduction's solver.
///
/// The transport at the end moves `round(k·(N'_a − N'_b))` per face as one
/// integer — conservation is independent of how well the solve converged.
pub fn diffuse_implicit(
    spec: GridSpec,
    field: &mut [i128],
    mobile: &[bool],
    k: f64,
    iters: u32,
    omega: f64,
    scratch: &mut DiffuseScratch,
) {
    diffuse_implicit_perm(spec, field, mobile, None, k, iters, omega, scratch)
}

/// Implicit diffusion with a per-cell permeability barrier.
///
/// Each face carries its own rate `k · p_a · p_b`, so the linear system stops
/// having a single coefficient and gains one per face. Diagonal dominance — and
/// therefore Gauss–Seidel's convergence at any step size — survives untouched,
/// because every off-diagonal entry is still non-negative and the diagonal is
/// still one plus their sum.
pub fn diffuse_implicit_perm(
    spec: GridSpec,
    field: &mut [i128],
    mobile: &[bool],
    perm: Option<&[f64]>,
    k: f64,
    iters: u32,
    omega: f64,
    scratch: &mut DiffuseScratch,
) {
    let cells = spec.cells();
    debug_assert_eq!(field.len(), cells);
    if cells == 0 || k <= 0.0 {
        return;
    }
    prepare(spec, field, mobile, scratch);

    for z in 0..spec.nz {
        for y in 0..spec.ny {
            for x in 0..spec.nx {
                let i = spec.idx(x, y, z);
                let v = field[i] as f64;
                scratch.n0[i] = v;
                scratch.n[i] = v;
                scratch.parity[i] = ((x + y + z) & 1) as u8;
            }
        }
    }

    // Neighbour lists are cheap to rebuild per call at these grid sizes; a
    // cached adjacency is an optimisation that waits for a profiler to ask.
    let neighbours = |x: usize, y: usize, z: usize| -> ([usize; 6], usize) {
        let mut nb = [0usize; 6];
        let mut c = 0;
        let i = spec.idx(x, y, z);
        let push = |j: usize, c: &mut usize, nb: &mut [usize; 6]| {
            if mobile[j] {
                nb[*c] = j;
                *c += 1;
            }
        };
        if mobile[i] {
            if x > 0 {
                push(spec.idx(x - 1, y, z), &mut c, &mut nb);
            }
            if x + 1 < spec.nx {
                push(spec.idx(x + 1, y, z), &mut c, &mut nb);
            }
            if y > 0 {
                push(spec.idx(x, y - 1, z), &mut c, &mut nb);
            }
            if y + 1 < spec.ny {
                push(spec.idx(x, y + 1, z), &mut c, &mut nb);
            }
            if z > 0 {
                push(spec.idx(x, y, z - 1), &mut c, &mut nb);
            }
            if z + 1 < spec.nz {
                push(spec.idx(x, y, z + 1), &mut c, &mut nb);
            }
        }
        (nb, c)
    };

    for _ in 0..iters {
        for colour in 0..2u8 {
            for z in 0..spec.nz {
                for y in 0..spec.ny {
                    for x in 0..spec.nx {
                        let i = spec.idx(x, y, z);
                        if scratch.parity[i] != colour || !mobile[i] {
                            continue;
                        }
                        let (nb, c) = neighbours(x, y, z);
                        if c == 0 {
                            continue;
                        }
                        let mut sum_kn = 0.0;
                        let mut sum_k = 0.0;
                        for &j in nb.iter().take(c) {
                            let kf = k * perm_of(perm, i) * perm_of(perm, j);
                            sum_kn += kf * scratch.n[j];
                            sum_k += kf;
                        }
                        let target = (scratch.n0[i] + sum_kn) / (1.0 + sum_k);
                        scratch.n[i] = (1.0 - omega) * scratch.n[i] + omega * target;
                    }
                }
            }
        }
    }

    // Residual, in counts: how far the final iterate is from actually solving
    // its own equation. Reported for the caller's judgement.
    let mut res: f64 = 0.0;
    for z in 0..spec.nz {
        for y in 0..spec.ny {
            for x in 0..spec.nx {
                let i = spec.idx(x, y, z);
                if !mobile[i] {
                    continue;
                }
                let (nb, c) = neighbours(x, y, z);
                if c == 0 {
                    continue;
                }
                let mut acc = 0.0;
                for &j in nb.iter().take(c) {
                    let kf = k * perm_of(perm, i) * perm_of(perm, j);
                    acc += kf * (scratch.n[j] - scratch.n[i]);
                }
                let r = (scratch.n[i] - scratch.n0[i]) - acc;
                if r.abs() > res {
                    res = r.abs();
                }
            }
        }
    }
    scratch.residual = res;

    // Transport: the only place matter moves, and it moves in integers —
    // each face's proposal capped by the donor's pocket.
    for_each_face(spec, |a, b| {
        if !mobile[a] || !mobile[b] {
            return;
        }
        let kf = k * perm_of(perm, a) * perm_of(perm, b);
        let raw = kf * (scratch.n[a] - scratch.n[b]);
        transfer(raw, a, b, scratch, field);
    });
}
