//! Vesicle geometry and inheritance without a genome — Phase 5i, step 1.
//!
//! Two anchors, both arithmetic rather than tuned:
//!
//!   * **v = 1/√2** is exactly the reduced volume at which a membrane's area is
//!     enough to close two equal spheres. Above it symmetric fission is
//!     geometrically impossible; below it, available.
//!   * **binomial variance**. Splitting `n` copies by a fair coin transmits them
//!     with a relative error of `1/(2√n)`. That is not a fidelity setting — it is
//!     what counting does — and it is why heredity gets more accurate as copy
//!     number rises.

use au_chem::vesicle::{
    composition_distance, fate, partition, reduced_volume, shape, Fate, FISSION_REDUCED_VOLUME,
};
use au_core::{Domain, Rng};

fn rng(k: u64) -> Rng {
    Rng::derive(90_210, Domain::CHEMISTRY, k, 0)
}

const A0: f64 = 4.0e-19; // headgroup area, m² — Phase 5f's constant

// ─── Geometry ────────────────────────────────────────────────────────────────

/// A sphere is the shape that encloses the most volume per unit area, and the
/// reduced volume is defined so that it scores exactly 1. If this drifts,
/// everything built on it is measuring something else.
#[test]
fn a_sphere_has_reduced_volume_one() {
    for r in [1.0e-9, 1.0e-6, 1.0e-3, 1.0] {
        let area = 4.0 * std::f64::consts::PI * r * r;
        let vol = (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
        let v = reduced_volume(area, vol).unwrap();
        assert!((v - 1.0).abs() < 1e-12, "r = {}: v = {}", r, v);
    }
}

/// **The anchor.** Two equal spheres holding the same total volume as one bag
/// need exactly √2 times its area — so the reduced volume of the two-sphere
/// shape is 1/√2, and that is the threshold, derived rather than chosen.
#[test]
fn two_equal_spheres_have_reduced_volume_one_over_root_two() {
    for r in [1.0e-9, 1.0e-6, 1.0] {
        let area = 2.0 * 4.0 * std::f64::consts::PI * r * r;
        let vol = 2.0 * (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
        let v = reduced_volume(area, vol).unwrap();
        assert!((v - FISSION_REDUCED_VOLUME).abs() < 1e-12, "r = {}: v = {}", r, v);
    }
    assert!((FISSION_REDUCED_VOLUME - 0.7071067811865476).abs() < 1e-15);
}

/// The reduced volume is a shape, not a size: scaling a bag up changes nothing
/// about what it can do. A threshold that moved with size would be a size
/// threshold wearing a disguise.
#[test]
fn the_criterion_is_scale_free() {
    let base = shape(1_000_000, 500_000, 1.0e24, A0).unwrap();
    for k in [2i128, 10, 1000] {
        // Area scales as k, so volume must scale as k^(3/2) to hold the shape.
        let v = shape(
            1_000_000 * k,
            (500_000.0 * (k as f64).powf(1.5)) as i128,
            1.0e24,
            A0,
        )
        .unwrap();
        assert!(
            (v.reduced_volume - base.reduced_volume).abs() / base.reduced_volume < 1e-6,
            "k = {}: {} vs {}",
            k,
            v.reduced_volume,
            base.reduced_volume
        );
    }
}

/// Osmotic volume is derived, not declared: dilute the surroundings and the bag
/// swells, because a membrane passes water and not solute.
#[test]
fn a_bag_swells_when_its_surroundings_are_diluted() {
    let concentrated = shape(1_000_000, 400_000, 2.0e24, A0).unwrap();
    let dilute = shape(1_000_000, 400_000, 1.0e24, A0).unwrap();
    assert!(dilute.volume_m3 > concentrated.volume_m3);
    assert!((dilute.volume_m3 / concentrated.volume_m3 - 2.0).abs() < 1e-9);
}

/// How many osmolytes a bag of `membrane` amphiphiles would hold if it were a
/// sphere at osmotic balance. Everything about a bag's fate is where it sits
/// relative to this number — which is itself derived, not chosen.
fn sphere_balance(membrane: i128, c_external: f64) -> i128 {
    let area = membrane as f64 * A0;
    let r = (area / (4.0 * std::f64::consts::PI)).sqrt();
    let v = (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
    (c_external * v) as i128
}

/// The trade-off, stated as a test.
///
/// A bag that makes osmolytes faster than membrane bursts; one that makes
/// membrane faster than osmolytes divides. **Neither fate mentions a size.**
/// What decides is where the bag sits relative to its own spherical balance —
/// the ratio of two chemistries, not a threshold anyone wrote down.
#[test]
fn fate_is_decided_by_a_ratio_of_two_chemistries_not_by_size() {
    let c_ext = 1.0e24;
    for membrane in [200_000i128, 1_800_000, 40_000_000] {
        let balance = sphere_balance(membrane, c_ext);

        // Surface to spare: half the solute a sphere would hold.
        let membranous = shape(membrane, balance / 2, c_ext, A0).unwrap();
        // Overfilled: twice what the area can enclose.
        let swollen = shape(membrane, balance * 2, c_ext, A0).unwrap();
        // Just under balance: a sphere, and it stays one.
        let poised = shape(membrane, balance * 90 / 100, c_ext, A0).unwrap();

        assert_eq!(fate(membranous.reduced_volume, 0.03), Fate::Fission, "N = {}", membrane);
        assert_eq!(fate(swollen.reduced_volume, 0.03), Fate::Lysis, "N = {}", membrane);
        assert_eq!(fate(poised.reduced_volume, 0.03), Fate::Intact, "N = {}", membrane);
    }
}

/// And the same ratio gives the same fate across two hundred–fold in size: a
/// bacterium-scale bag and a speck are governed by one rule.
#[test]
fn two_bags_of_very_different_sizes_share_a_fate_if_they_share_a_ratio() {
    let c_ext = 1.0e24;
    let small = shape(200_000, sphere_balance(200_000, c_ext) / 2, c_ext, A0).unwrap();
    let large = shape(40_000_000, sphere_balance(40_000_000, c_ext) / 2, c_ext, A0).unwrap();
    assert_eq!(fate(small.reduced_volume, 0.03), fate(large.reduced_volume, 0.03));
    assert!((small.reduced_volume - large.reduced_volume).abs() < 1e-3);
}

/// A bilayer tolerates a few percent of areal strain before it fails, so the
/// lysis line sits slightly above the sphere rather than exactly on it.
#[test]
fn a_membrane_stretches_a_little_before_it_bursts() {
    assert_eq!(fate(1.01, 0.03), Fate::Intact, "3% strain should survive 1% overfill");
    assert_eq!(fate(1.20, 0.03), Fate::Lysis);
    assert_eq!(fate(1.01, 0.0), Fate::Lysis, "with no tolerance, v > 1 must fail");
}

/// Right at the line, fission is available — the bound is inclusive, because at
/// exactly 1/√2 the two spheres close exactly.
#[test]
fn the_fission_bound_is_inclusive() {
    assert_eq!(fate(FISSION_REDUCED_VOLUME, 0.03), Fate::Fission);
    assert_eq!(fate(FISSION_REDUCED_VOLUME + 1e-9, 0.03), Fate::Intact);
}

// ─── Partition ───────────────────────────────────────────────────────────────

/// Division moves molecules. It does not make or destroy them, and a protocell
/// that could would be a free-matter machine of exactly the kind Phase 5h's
/// ledger exists to catch.
#[test]
fn division_never_costs_a_molecule() {
    let mut r = rng(1);
    let parent: Vec<i128> = vec![1, 2, 7, 99, 1_000, 65_536, 5_000_000, 0, 3];
    for _ in 0..200 {
        let (a, b) = partition(&parent, 0.5, &mut r);
        for i in 0..parent.len() {
            assert_eq!(a[i] + b[i], parent[i], "species {} lost or gained", i);
            assert!(a[i] >= 0 && b[i] >= 0);
        }
    }
}

/// Same seed, same split — always. Stochastic is not the same thing as
/// nondeterministic, and this project's whole foundation depends on the
/// difference.
#[test]
fn the_same_stream_gives_the_same_daughters() {
    let parent: Vec<i128> = vec![10, 500, 20_000, 3_000_000];
    let (a1, b1) = partition(&parent, 0.5, &mut rng(7));
    let (a2, b2) = partition(&parent, 0.5, &mut rng(7));
    assert_eq!(a1, a2);
    assert_eq!(b1, b2);
    let (a3, _) = partition(&parent, 0.5, &mut rng(8));
    assert_ne!(a1, a3, "two different streams produced identical draws");
}

/// **The quantitative anchor.** A fair split of `n` copies has standard
/// deviation `√n/2`. Measured over many divisions, not asserted.
#[test]
fn the_split_has_binomial_variance() {
    let mut r = rng(11);
    for n in [64i128, 1_000, 250_000] {
        let trials = 4_000;
        let (mut sum, mut sumsq) = (0.0f64, 0.0f64);
        for _ in 0..trials {
            let (a, _) = partition(&[n], 0.5, &mut r);
            sum += a[0] as f64;
            sumsq += (a[0] as f64).powi(2);
        }
        let mean = sum / trials as f64;
        let sd = (sumsq / trials as f64 - mean * mean).sqrt();
        let expect_sd = (n as f64).sqrt() / 2.0;
        assert!(
            (mean - n as f64 / 2.0).abs() < 4.0 * expect_sd / (trials as f64).sqrt(),
            "n = {}: mean {} is not n/2",
            n,
            mean
        );
        assert!(
            (sd / expect_sd - 1.0).abs() < 0.08,
            "n = {}: sd {:.2}, expected {:.2}",
            n,
            sd,
            expect_sd
        );
    }
}

/// An uneven pinch gives an uneven split, in proportion — so a bag that divides
/// asymmetrically passes on its composition asymmetrically too.
#[test]
fn an_uneven_pinch_splits_in_proportion() {
    let mut r = rng(13);
    let (a, b) = partition(&[10_000_000], 0.25, &mut r);
    let frac = a[0] as f64 / (a[0] + b[0]) as f64;
    assert!((frac - 0.25).abs() < 0.01, "asked for 0.25, got {:.4}", frac);
}

// ─── Inheritance without a genome ────────────────────────────────────────────

/// **The result this step exists to demonstrate.**
///
/// There is no genome here. No template, no copying, no sequence — only a bag
/// splitting in two. And yet a daughter resembles its parent far more closely
/// than it resembles an arbitrary composition of the same size, because the
/// coin is fair *per molecule* and the proportions therefore survive the split.
///
/// This is compositional inheritance (Segré, Ben-Eli & Lancet, 2000): heredity
/// that is carried by *what a thing is made of* rather than by anything written
/// down. It is the weakest possible form of the property, and it is real.
#[test]
fn a_daughter_resembles_its_parent_without_anything_being_copied() {
    let mut r = rng(21);
    // A composition with structure: some species abundant, some rare.
    let parent: Vec<i128> =
        (0..40).map(|i| 1_000_000 / (1 + i as i128 * i as i128).max(1)).collect();

    let (daughter, _) = partition(&parent, 0.5, &mut r);
    let inherited = composition_distance(&parent, &daughter);

    // A stranger of the same total size, built by shuffling the same abundances
    // onto different species.
    let mut stranger = parent.clone();
    stranger.reverse();
    let unrelated = composition_distance(&parent, &stranger);

    assert!(
        inherited < unrelated / 50.0,
        "daughter distance {:.6} vs stranger {:.6} — nothing was inherited",
        inherited,
        unrelated
    );
}

/// **Fidelity is copy number.** The same split is far more faithful when the
/// molecules are abundant, because the error is `1/(2√n)`. Nothing sets this;
/// it is what counting does — and it is the reason a real cell keeps one genome
/// rather than one molecule of every protein.
#[test]
fn heredity_gets_more_faithful_as_copy_number_rises() {
    let mut r = rng(31);
    let mut previous = f64::INFINITY;
    for scale in [10i128, 1_000, 100_000, 10_000_000] {
        // Same recipe at every scale: only the absolute counts change.
        let parent: Vec<i128> = (1..=20).map(|i| scale * i as i128).collect();
        // Average over several divisions so the comparison is of typical error,
        // not of one lucky draw.
        let mut total = 0.0;
        let trials = 200;
        for _ in 0..trials {
            let (d, _) = partition(&parent, 0.5, &mut r);
            total += composition_distance(&parent, &d);
        }
        let err = total / trials as f64;
        assert!(err < previous, "fidelity did not improve at scale {}: {:.3e}", scale, err);
        previous = err;
    }
    // And the best case is very faithful indeed.
    assert!(previous < 1.0e-8, "even abundant species were transmitted badly: {:.3e}", previous);
}

/// The honest limit, written as a test so it cannot be quietly forgotten.
///
/// Compositional inheritance decays. Take a daughter, let it grow back to its
/// parent's size, divide again, and repeat: the distance from the founder grows
/// with every generation, because each split adds its own √n noise and nothing
/// corrects it. There is no proofreading, no template to compare against — a
/// composition cannot be *checked*, only re-drawn.
///
/// This is the known weakness of compositional genomes (Vasas, Szathmáry &
/// Santos, 2010): they carry information but cannot hold it still, so selection
/// on them has very little to grip. Recording it here means the next phase knows
/// what problem it is solving rather than discovering it later.
#[test]
fn compositional_inheritance_decays_across_generations() {
    let mut r = rng(41);
    let founder: Vec<i128> = (1..=25).map(|i| 40_000 * i as i128).collect();
    let mut lineage = founder.clone();
    let mut drift = Vec::new();
    for _ in 0..12 {
        let (d, _) = partition(&lineage, 0.5, &mut r);
        // Grow back: the daughter doubles what it has, which is exactly the
        // point — growth amplifies whatever the split happened to give it.
        lineage = d.iter().map(|&x| x * 2).collect();
        drift.push(composition_distance(&founder, &lineage));
    }
    assert!(
        drift[11] > drift[0],
        "no drift at all across 12 generations: {:.3e} -> {:.3e}",
        drift[0],
        drift[11]
    );
}
