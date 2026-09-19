//! Collision-theory prefactors, checked against the numbers kinetics textbooks
//! quote. These are the *anchors* for Phase 5e: not internal consistency, but
//! agreement with measured chemistry.

use au_chem::kinetics::*;


/// **The Eyring frequency.** A unimolecular attempt frequency is `k_B·T/h`,
/// which is ~2×10¹³ s⁻¹ at 1000 K — the figure quoted for bond vibration, and
/// the reason unimolecular A-factors are always reported near 10¹³.
#[test]
fn a_unimolecular_prefactor_is_a_vibration_frequency() {
    let r = KineticRules::default();
    let a = prefactor_at(unimolecular_prefactor(&r), 1, 1000.0);
    assert!(
        (1.0e13..4.0e13).contains(&a),
        "unimolecular A at 1000 K was {:.3e} s^-1, expected ~2e13",
        a
    );
    // Linear in temperature, exactly.
    let hot = prefactor_at(unimolecular_prefactor(&r), 1, 2000.0);
    assert!((hot / a - 2.0).abs() < 1.0e-12, "kT/h must be linear in T");
}

/// **The gas-kinetic collision rate.** Two small molecules at room temperature
/// collide at ~10⁻¹⁶ m³·molecule⁻¹·s⁻¹, which is ~10¹¹ L·mol⁻¹·s⁻¹ — the
/// canonical bimolecular A-factor, and the ceiling that diffusion-limited
/// reactions approach.
#[test]
fn a_bimolecular_prefactor_is_the_gas_kinetic_collision_rate() {
    let r = KineticRules::default();
    // Two ~30 Da species of three atoms each.
    let structural = bimolecular_prefactor(30_000, 3, 30_000, 3, &r);
    let a = prefactor_at(structural, 2, 300.0);
    assert!(
        (1.0e-17..1.0e-15).contains(&a),
        "bimolecular A at 300 K was {:.3e} m^3/molecule/s, expected ~1e-16",
        a
    );
    // Expressed the way a chemist would read it.
    let per_mol_per_litre = a * 6.022_140_76e23 * 1000.0;
    assert!(
        (1.0e10..1.0e12).contains(&per_mol_per_litre),
        "that is {:.3e} L/mol/s, expected ~1e11",
        per_mol_per_litre
    );
}

/// Light molecules collide more often: the rate carries `1/√μ`, so hydrogen is
/// far quicker off the mark than anything heavy. Same law that set diffusion
/// rates in Phase 5b, arriving here independently.
#[test]
fn lighter_molecules_collide_more_often() {
    let r = KineticRules::default();
    let light = bimolecular_prefactor(2_016, 2, 2_016, 2, &r);
    let heavy = bimolecular_prefactor(32_000, 2, 32_000, 2, &r);
    let ratio = light / heavy;
    let expected = (32_000.0f64 / 2_016.0).sqrt(); // √(μ_heavy/μ_light)
    assert!(
        (ratio / expected - 1.0).abs() < 0.02,
        "collision ratio {:.3} vs √mass ratio {:.3}",
        ratio,
        expected
    );
}

/// Bigger molecules present a bigger target: the cross-section grows as the
/// square of the summed radii, and radius as the cube root of atom count.
#[test]
fn bigger_molecules_present_a_bigger_target() {
    let r = KineticRules::default();
    // Same masses, different sizes, so only σ differs.
    let small = bimolecular_prefactor(30_000, 1, 30_000, 1, &r);
    let big = bimolecular_prefactor(30_000, 8, 30_000, 8, &r);
    let ratio = big / small;
    let expected = ((8f64.cbrt() + 8f64.cbrt()) / 2.0).powi(2);
    assert!(
        (ratio / expected - 1.0).abs() < 0.02,
        "cross-section ratio {:.3} vs geometric {:.3}",
        ratio,
        expected
    );
}

/// Collision rates rise as √T — weakly, which is why the exponential dominates
/// and why Arrhenius plots are straight enough to be useful.
#[test]
fn collision_rates_rise_as_the_square_root_of_temperature() {
    let r = KineticRules::default();
    let s = bimolecular_prefactor(20_000, 3, 20_000, 3, &r);
    let cold = prefactor_at(s, 2, 300.0);
    let hot = prefactor_at(s, 2, 1200.0);
    assert!((hot / cold - 2.0).abs() < 1.0e-9, "a fourfold T should double the rate");
}

/// **Three-body encounters are rare**, by the volume factor. This is the
/// physical reason behind Phase 5c's empirical finding that a catalyst only pays
/// its way when the barrier it removes is large: it must find a third body, and
/// that costs orders of magnitude before any barrier is considered.
#[test]
fn termolecular_encounters_are_orders_of_magnitude_rarer() {
    let r = KineticRules::default();
    let two = prefactor_at(bimolecular_prefactor(30_000, 3, 30_000, 3, &r), 2, 1000.0);
    let three = prefactor_at(termolecular_prefactor(30_000, 3, 30_000, 3, &r), 3, 1000.0);
    let penalty = two / three;
    assert!(
        penalty > 1.0e20,
        "a third body should cost many orders of magnitude, got {:.2e}",
        penalty
    );
    assert!((penalty - 1.0 / r.encounter_volume_m3).abs() / penalty < 1.0e-9);
}

/// The steric factor is the only dial, and it does exactly what it says:
/// scales every prefactor, changing rates and never equilibria.
#[test]
fn the_steric_factor_is_the_only_free_parameter() {
    let strict = KineticRules { steric_factor: 0.01, ..Default::default() };
    let loose = KineticRules::default();
    let a = bimolecular_prefactor(30_000, 3, 30_000, 3, &strict);
    let b = bimolecular_prefactor(30_000, 3, 30_000, 3, &loose);
    assert!((a / b - 0.01).abs() < 1.0e-12);
}

/// Prefactors are a property of the molecules, so the same pair always gives the
/// same answer whichever way round it is asked.
#[test]
fn collision_prefactors_are_symmetric_and_reproducible() {
    let r = KineticRules::default();
    let ab = bimolecular_prefactor(18_015, 3, 44_010, 3, &r);
    let ba = bimolecular_prefactor(44_010, 3, 18_015, 3, &r);
    assert_eq!(ab, ba, "a collision is not symmetric in its partners");
    assert_eq!(ab, bimolecular_prefactor(18_015, 3, 44_010, 3, &r));
}
