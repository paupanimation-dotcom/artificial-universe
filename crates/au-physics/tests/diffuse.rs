//! Diffusion, held to the house standard: exact conservation regardless of
//! solver quality, agreement between schemes where both are valid, stability
//! past the explicit cliff, and one analytic anchor — the Gaussian.
//!
//! The anchor: a point release diffusing in one dimension spreads as a Gaussian
//! whose variance grows linearly, `σ² = 2·D·t` (Einstein, 1905 — the Brownian
//! motion paper). In lattice units with rate `k = D·dt/dx²`, that is
//! `σ² = 2·k·steps` measured in cells². The engine is not told this; it is told
//! only "flux is proportional to the difference across a face", and the
//! Gaussian is what that rule *does*.

use au_physics::*;

fn line(n: usize) -> GridSpec {
    GridSpec::new(n, 1, 1, 1.0)
}

fn all_mobile(n: usize) -> Vec<bool> {
    vec![true; n]
}

fn total(f: &[i128]) -> i128 {
    f.iter().sum()
}

// ═══ Conservation and basic sanity ══════════════════════════════════════════

/// The sum is invariant to the last count, whatever the rate and pattern.
#[test]
fn explicit_diffusion_conserves_exactly() {
    let spec = line(41);
    let mut f = vec![0i128; 41];
    f[7] = 1_000_003;
    f[20] = 999_999_999_999;
    f[33] = 17;
    let before = total(&f);
    let mobile = all_mobile(41);
    let mut scratch = DiffuseScratch::new();
    for _ in 0..500 {
        diffuse_explicit(spec, &mut f, &mobile, 0.2, &mut scratch);
    }
    assert_eq!(total(&f), before, "explicit diffusion leaked counts");
    assert!(f.iter().all(|&v| v >= 0), "a population went negative");
}

/// A sealed tube equilibrates — down to the quantisation floor, which is the
/// honest fixed point. With k = 0.25, a face whose difference is 1 rounds its
/// flux to zero, so the terminal state is not "every cell equal" but "no two
/// neighbours differ by more than 1" — a molecule-level staircase the module
/// docs promise and this test pins as the spec.
#[test]
fn diffusion_equilibrates_to_the_quantisation_floor() {
    let spec = line(16);
    let mut f = vec![0i128; 16];
    f[0] = 1_600_000;
    let mobile = all_mobile(16);
    let mut scratch = DiffuseScratch::new();
    for _ in 0..20_000 {
        diffuse_explicit(spec, &mut f, &mobile, 0.25, &mut scratch);
    }
    for i in 0..15 {
        assert!(
            (f[i] - f[i + 1]).abs() <= 1,
            "neighbours {} and {} still differ by {} — above the deadband",
            i,
            i + 1,
            f[i] - f[i + 1]
        );
    }
    let mean = total(&f) / 16;
    assert!(
        f.iter().all(|&v| (v - mean).abs() <= 16),
        "the staircase drifted further from the mean than the tube is long"
    );
}

/// Counts flow down their own gradient: the loaded cell loses, its empty
/// neighbour gains, in a single step.
#[test]
fn counts_flow_downhill_only() {
    let spec = line(3);
    let mut f = vec![0, 1_000_000, 0];
    let mobile = all_mobile(3);
    let mut scratch = DiffuseScratch::new();
    diffuse_explicit(spec, &mut f, &mobile, 0.2, &mut scratch);
    assert!(f[1] < 1_000_000, "the peak did not lose");
    assert!(f[0] > 0 && f[2] > 0, "the neighbours did not gain");
    assert_eq!(f[0], f[2], "a symmetric start must stay symmetric");
}

// ═══ The analytic anchor ════════════════════════════════════════════════════

/// A point release spreads as a Gaussian: variance grows as `2·k` per step
/// (Einstein, 1905). Measured over a long run, the engine's variance must track
/// the line — this is diffusion *being* diffusion, not merely conserving.
#[test]
fn a_point_release_spreads_with_einstein_variance() {
    let spec = line(201);
    let n = 201usize;
    let mut f = vec![0i128; n];
    f[100] = 1_000_000_000_000; // large, so integer rounding is invisible
    let mobile = all_mobile(n);
    let mut scratch = DiffuseScratch::new();
    let k = 0.2;
    let steps = 2000u32;
    for _ in 0..steps {
        diffuse_explicit(spec, &mut f, &mobile, k, &mut scratch);
    }
    // Variance about the centre, in cells².
    let tot = total(&f) as f64;
    let mut mean = 0.0;
    for (i, &v) in f.iter().enumerate() {
        mean += i as f64 * v as f64;
    }
    mean /= tot;
    let mut var = 0.0;
    for (i, &v) in f.iter().enumerate() {
        let d = i as f64 - mean;
        var += d * d * v as f64;
    }
    var /= tot;

    let predicted = 2.0 * k * steps as f64;
    let ratio = var / predicted;
    assert!(
        (0.97..1.03).contains(&ratio),
        "variance {:.1} vs Einstein's 2kt = {:.1} (ratio {:.4})",
        var,
        predicted,
        ratio
    );
}

// ═══ The implicit twin ══════════════════════════════════════════════════════

/// Where both schemes are valid (small k), they agree — same rule, same
/// face couplings, so any daylight between them means one is wrong.
#[test]
fn implicit_agrees_with_explicit_at_small_rates() {
    let n = 30;
    let spec = line(n);
    let mobile = all_mobile(n);
    let start = |f: &mut Vec<i128>| {
        f.clear();
        f.resize(n, 0);
        f[5] = 40_000_000;
        f[22] = 90_000_000;
    };
    let mut a = Vec::new();
    start(&mut a);
    let mut b = a.clone();
    let mut scratch = DiffuseScratch::new();
    for _ in 0..300 {
        diffuse_explicit(spec, &mut a, &mobile, 0.05, &mut scratch);
        diffuse_implicit(spec, &mut b, &mobile, 0.05, 60, 1.0, &mut scratch);
    }
    let tot = total(&a) as f64;
    for i in 0..n {
        let diff = (a[i] - b[i]).abs() as f64;
        assert!(
            diff / tot < 2.0e-3,
            "cell {}: explicit {} vs implicit {} — schemes disagree where both are valid",
            i,
            a[i],
            b[i]
        );
    }
}

/// A rate a thousand times past the explicit ceiling: bounded, sane, inside
/// the initial range, and arriving at the right destination. Explicit at this
/// k tears itself apart in a step; implicit shrugs. This is deep time for
/// chemistry.
///
/// One honest cost is visible here: the donor's-pocket limiter caps each face
/// at the donor's holdings per step, so an extreme stride spreads matter as a
/// front — about one cell per face per step — rather than teleporting to the
/// solved profile. Guaranteed non-negativity is bought with rate-limited
/// arrival; the run is long enough for the front to cross the tube and settle.
#[test]
fn implicit_is_stable_far_past_the_explicit_limit() {
    let n = 30;
    let spec = line(n);
    let mobile = all_mobile(n);
    let mut f = vec![0i128; n];
    f[10] = 5_000_000;
    let (lo, hi) = (0i128, 5_000_000i128);
    let mut scratch = DiffuseScratch::new();
    for _ in 0..200 {
        diffuse_implicit(spec, &mut f, &mobile, 250.0, 60, 1.0, &mut scratch);
    }
    for (i, &v) in f.iter().enumerate() {
        assert!(
            v >= lo && v <= hi,
            "cell {} left the initial range at {} — maximum-principle violation",
            i,
            v
        );
    }
    // And it should be nearly uniform: k=250 is far past mixing time.
    let mean = total(&f) / n as i128;
    assert!(
        f.iter().all(|&v| (v - mean).abs() < mean / 5),
        "at k=250 the tube should be essentially mixed"
    );
}

/// Conservation is independent of how well the solve converged: one sweep, an
/// enormous rate, and the books still balance to the last count. Wrong physics
/// is permitted; lost matter is not.
#[test]
fn conservation_does_not_depend_on_convergence() {
    let n = 25;
    let spec = line(n);
    let mobile = all_mobile(n);
    for iters in [1u32, 2, 5, 60] {
        let mut f = vec![0i128; n];
        f[3] = 123_456_789;
        f[19] = 987_654_321;
        let before = total(&f);
        let mut scratch = DiffuseScratch::new();
        for _ in 0..100 {
            diffuse_implicit(spec, &mut f, &mobile, 300.0, iters, 1.0, &mut scratch);
        }
        assert_eq!(
            total(&f),
            before,
            "a {}-sweep solve leaked {} counts",
            iters,
            total(&f) - before
        );
        assert!(f.iter().all(|&v| v >= 0), "{} sweeps drove a count negative", iters);
    }
}

// ═══ The mobility mask ══════════════════════════════════════════════════════

/// Solids and vacuum are walls: an immobile cell neither gives nor takes, and
/// the mobile cells conserve among themselves. This is the rule that lets a
/// freezing ocean trap its chemistry in place.
#[test]
fn immobile_cells_are_walls() {
    let n = 9;
    let spec = line(n);
    let mut mobile = all_mobile(n);
    mobile[4] = false; // a frozen cell splits the tube in two
    let mut f = vec![0i128; n];
    f[1] = 800_000; // left basin
    f[4] = 55_555; // trapped in the ice
    f[7] = 200_000; // right basin
    let left_before: i128 = f[..4].iter().sum();
    let right_before: i128 = f[5..].iter().sum();

    let mut scratch = DiffuseScratch::new();
    for _ in 0..2000 {
        diffuse_explicit(spec, &mut f, &mobile, 0.25, &mut scratch);
        diffuse_implicit(spec, &mut f, &mobile, 3.0, 40, 1.0, &mut scratch);
    }

    assert_eq!(f[4], 55_555, "the frozen cell's cargo moved");
    assert_eq!(f[..4].iter().sum::<i128>(), left_before, "the left basin leaked past the ice");
    assert_eq!(f[5..].iter().sum::<i128>(), right_before, "the right basin leaked past the ice");
    // Each basin equilibrated internally.
    assert!((f[0] - f[3]).abs() <= 2, "left basin failed to mix");
    assert!((f[5] - f[8]).abs() <= 2, "right basin failed to mix");
}

/// Same input, same output, to the bit — fixed sweeps in a fixed order, no
/// state-dependent work.
#[test]
fn diffusion_is_bit_reproducible() {
    let n = 33;
    let spec = line(n);
    let mobile = all_mobile(n);
    let run = || {
        let mut f = vec![0i128; n];
        f[8] = 77_000_000;
        f[25] = 13_000_000;
        let mut scratch = DiffuseScratch::new();
        for _ in 0..200 {
            diffuse_explicit(spec, &mut f, &mobile, 0.15, &mut scratch);
            diffuse_implicit(spec, &mut f, &mobile, 4.0, 30, 1.0, &mut scratch);
        }
        f
    };
    assert_eq!(run(), run(), "diffusion is not reproducible");
}

/// The deep-time regression, locked: at k = 1000 on a 41-cell tube, the sweep
/// budget from `sweeps_for_deep_time` (2·n²) converges the solve so completely
/// that after the pocket-limited front has crossed and settled, the tube sits
/// at the *exact* uniform count — worst deviation zero. The first attempt used
/// 24 sweeps here; the k-amplified solver error manufactured a checkerboard,
/// starved the end cells to literal zero, and taught the module its "How many
/// sweeps" section. This test keeps that lesson executable.
#[test]
fn deep_time_sweeps_scale_with_span_squared_and_mix_exactly() {
    let n = 41;
    let spec = line(n);
    let mobile = all_mobile(n);
    let mut f = vec![0i128; n];
    f[20] = 41_000_000;
    let sweeps = sweeps_for_deep_time(spec, 24);
    assert_eq!(sweeps, 2 * 41 * 41, "the sweep law changed silently");
    let mut scratch = DiffuseScratch::new();
    for _ in 0..60 {
        diffuse_implicit(spec, &mut f, &mobile, 1000.0, sweeps, 1.0, &mut scratch);
    }
    assert_eq!(total(&f), 41_000_000, "deep time leaked");
    for (i, &v) in f.iter().enumerate() {
        assert_eq!(v, 1_000_000, "cell {} at {} — not the exact uniform count", i, v);
    }
    assert!(scratch.residual < 1.0, "the solve should be essentially exact here");
}

// ═══ Permeability barriers ══════════════════════════════════════════════════

/// A sealed cell keeps what it has and takes nothing in. Permeability zero is
/// the limit case of a membrane, and it must be exact: not "almost nothing
/// crosses" but *nothing*, so a fully closed compartment is genuinely closed.
#[test]
fn a_sealed_cell_neither_gives_nor_takes() {
    let n = 9;
    let spec = line(n);
    let mobile = all_mobile(n);
    let mut perm = vec![1.0; n];
    perm[4] = 0.0;
    let mut f = vec![0i128; n];
    f[1] = 800_000;
    f[4] = 55_555;
    f[7] = 200_000;
    let mut scratch = DiffuseScratch::new();
    for _ in 0..2000 {
        diffuse_explicit_perm(spec, &mut f, &mobile, Some(&perm), 0.25, &mut scratch);
        diffuse_implicit_perm(spec, &mut f, &mobile, Some(&perm), 3.0, 40, 1.0, &mut scratch);
    }
    assert_eq!(f[4], 55_555, "a sealed cell exchanged contents");
    assert_eq!(f[..4].iter().sum::<i128>(), 800_000, "matter crossed the seal");
    assert_eq!(f[5..].iter().sum::<i128>(), 200_000, "matter crossed the seal");
}

/// A partial barrier slows exchange without stopping it — which is what a real
/// membrane does, and why a compartment is a cell rather than a tomb. Same tube,
/// same time, three permeabilities: the more closed the cell, the more it keeps.
#[test]
fn a_barrier_slows_exchange_in_proportion_to_its_permeability() {
    let n = 11;
    let spec = line(n);
    let mobile = all_mobile(n);
    let retained = |p: f64| -> i128 {
        let mut perm = vec![1.0; n];
        perm[5] = p;
        let mut f = vec![0i128; n];
        f[5] = 10_000_000;
        let mut scratch = DiffuseScratch::new();
        for _ in 0..400 {
            diffuse_explicit_perm(spec, &mut f, &mobile, Some(&perm), 0.2, &mut scratch);
        }
        f[5]
    };
    let (open, half, tight) = (retained(1.0), retained(0.3), retained(0.02));
    assert!(tight > half && half > open, "retention {} / {} / {}", tight, half, open);
    assert!(open > 0);
    assert!(tight < 10_000_000, "a partial barrier must leak eventually");
}

/// **A barrier is a membrane, not a pump.** It resists a molecule leaving
/// exactly as it resists one arriving; if it did not, a world full of membranes
/// would fill or empty itself for free.
///
/// Tested as an exact mirror rather than a ratio. Load each end of a barriered
/// tube in turn: if the operator is symmetric, the second profile must be the
/// first one reversed, to the molecule. (An earlier version compared outward
/// against inward flux from unequal starting gradients, measured a ratio of
/// exactly 6.0, and the gradients differed by exactly 6 — it was measuring its
/// own setup.)
#[test]
fn a_barrier_is_a_membrane_not_a_pump() {
    let n = 5;
    let spec = line(n);
    let mobile = all_mobile(n);
    let perm = vec![1.0, 1.0, 0.15, 1.0, 1.0];
    let run = |load_at: usize, p: &[f64]| -> Vec<i128> {
        let mut f = vec![0i128; n];
        f[load_at] = 1_000_000;
        let mut scratch = DiffuseScratch::new();
        for _ in 0..60 {
            diffuse_explicit_perm(spec, &mut f, &mobile, Some(p), 0.2, &mut scratch);
        }
        f
    };
    let left = run(0, &perm);
    let mut right = run(4, &perm);
    right.reverse();
    assert_eq!(left, right, "the barrier treats the two directions differently");

    let open = run(0, &vec![1.0; n]);
    assert!(
        left[4] < open[4] / 2,
        "barrier {} vs open {} at the far end — nothing was impeded",
        left[4],
        open[4]
    );
}

/// Permeability changes rates, never totals. Whatever the barrier profile and
/// whichever solver, the books balance to the last molecule.
#[test]
fn permeability_never_costs_a_molecule() {
    let n = 13;
    let spec = line(n);
    let mobile = all_mobile(n);
    let perm: Vec<f64> = (0..n).map(|i| if i % 3 == 0 { 0.05 } else { 1.0 }).collect();
    let mut f = vec![0i128; n];
    f[2] = 987_654_321;
    f[9] = 123_456_789;
    let before: i128 = f.iter().sum();
    let mut scratch = DiffuseScratch::new();
    for _ in 0..500 {
        diffuse_explicit_perm(spec, &mut f, &mobile, Some(&perm), 0.15, &mut scratch);
        diffuse_implicit_perm(spec, &mut f, &mobile, Some(&perm), 200.0, 400, 1.0, &mut scratch);
    }
    assert_eq!(f.iter().sum::<i128>(), before, "a barrier leaked matter");
    assert!(f.iter().all(|&v| v >= 0), "a barrier drove a population negative");
}
