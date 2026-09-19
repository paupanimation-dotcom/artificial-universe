//! Where reaction rates actually come from.
//!
//! # The correction this module exists to make
//!
//! Until now every generated reaction carried the *same* declared
//! pre-exponential factor — one number, applied to unimolecular and bimolecular
//! reactions alike. That cannot be right, and not merely as a matter of taste:
//! the two have **different units**. Reading the reactor's own rate law,
//!
//! ```text
//!     events = k · Π nᵢ · V · dt
//! ```
//!
//! a unimolecular `k` must be s⁻¹ (an attempt frequency) while a bimolecular one
//! must be m³·molecule⁻¹·s⁻¹ (a swept volume per second). A single number cannot
//! be both, so the old constant was not a physical quantity at all — it was a
//! fitting parameter that happened to make chemistry proceed at the populations
//! the demos used. It worked because a second error cancelled it: the engine ran
//! at concentrations far below anything real, and a prefactor far above anything
//! real brought the product back into range.
//!
//! Two compensating errors are fine for as long as nothing else has to agree
//! with either of them. That ended the moment membranes arrived carrying genuine
//! numbers — a millimolar critical concentration, a 0.4 nm² headgroup — because
//! a calibrated layer cannot sit on an uncalibrated one. So the prefactor is
//! derived here instead, from collision theory, out of quantities the molecular
//! graphs already carry.
//!
//! # What is derived, and from what
//!
//! **Unimolecular** (a molecule falling apart, or a bond order shifting): the
//! attempt frequency is the Eyring factor `k_B·T/h` — how often a vibration
//! tries the barrier at all. About 2×10¹³ s⁻¹ at 1000 K, which is the textbook
//! figure. Nothing structural enters; only a dimensionless steric factor.
//!
//! **Bimolecular** (two molecules meeting): the collision frequency
//!
//! ```text
//!     A(T) = P · σ · √(8·k_B·T / π·μ)
//! ```
//!
//! where `σ = π(r_A + r_B)²` is the collision cross-section, `μ` the reduced
//! mass, and `P` a steric factor for the fraction of collisions correctly
//! oriented. Both `σ` and `μ` come off the graphs: mass from the summed atomic
//! masses, radius from the atom count as `r₀·n^(1/3)` — a molecule of `n` atoms
//! occupies roughly `n` atomic volumes, so its radius grows as the cube root.
//! This lands at ~10⁻¹⁶ m³·molecule⁻¹·s⁻¹, i.e. ~10¹¹ L·mol⁻¹·s⁻¹, which is the
//! gas-kinetic value every kinetics textbook quotes.
//!
//! **Termolecular** (the catalysed reactions of Phase 5c): a bimolecular
//! encounter that must find a third body inside a small encounter volume. The
//! extra factor of volume is why three-body reactions are rare, and why a
//! catalyst has to buy a large barrier reduction to be worth the trouble — which
//! is exactly what Phase 5c discovered empirically before there was any physical
//! reason for it in the code.
//!
//! # The consequence nobody has to pay for
//!
//! Forward and reverse reactions now have *different* prefactors, because an
//! association is bimolecular forwards and unimolecular backwards. Their ratio
//! is not noise — it is the reaction entropy:
//!
//! ```text
//!     K = (A_f / A_r) · exp(−ΔH/RT)
//! ```
//!
//! and `A_f/A_r` for an association is a very small volume, which says that
//! joining two molecules into one costs translational entropy. That term was
//! listed as a known omission for three phases. It is not being added here; it
//! *falls out* of counting molecules on each side, once the prefactors are
//! allowed to mean what they physically mean.

/// Boltzmann's constant, J·K⁻¹.
pub const BOLTZMANN: f64 = 1.380_649e-23;
/// Planck's constant, J·s.
pub const PLANCK: f64 = 6.626_070_15e-34;
/// One dalton, kg.
pub const DALTON: f64 = 1.660_539_066_60e-27;

/// The declared scales of collision kinetics — all dimensionless or geometric,
/// none of them a fitted rate.
#[derive(Clone, Copy, Debug)]
pub struct KineticRules {
    /// Fraction of collisions with the right orientation to react. Genuinely a
    /// fudge in collision theory too, which is why it is the one dial here:
    /// everything else is measured off the molecule.
    pub steric_factor: f64,
    /// Radius of a single atom, m. A molecule of `n` atoms is taken to have
    /// radius `r₀·n^(1/3)`.
    pub atom_radius_m: f64,
    /// The volume within which a third body counts as present for a
    /// termolecular encounter, m³.
    pub encounter_volume_m3: f64,
}

impl Default for KineticRules {
    fn default() -> Self {
        KineticRules {
            steric_factor: 1.0,
            atom_radius_m: 1.8e-10,
            encounter_volume_m3: 1.0e-28,
        }
    }
}

/// Molecular mass in kg, from millidaltons.
#[inline]
pub fn mass_kg(mass_mda: i64) -> f64 {
    mass_mda as f64 * 1.0e-3 * DALTON
}

/// Collision radius of a molecule with `n` atoms, m.
#[inline]
pub fn radius_m(atom_count: usize, rules: &KineticRules) -> f64 {
    rules.atom_radius_m * (atom_count.max(1) as f64).cbrt()
}

/// The structural part of a **unimolecular** prefactor.
///
/// Returns the dimensionless steric factor; the temperature-dependent Eyring
/// frequency `k_B·T/h` is applied at rate time, since it is not a property of
/// the molecule.
#[inline]
pub fn unimolecular_prefactor(rules: &KineticRules) -> f64 {
    rules.steric_factor
}

/// The structural part of a **bimolecular** prefactor, such that
/// `A(T) = returned · √T` in m³·molecule⁻¹·s⁻¹.
pub fn bimolecular_prefactor(
    mass_a_mda: i64,
    atoms_a: usize,
    mass_b_mda: i64,
    atoms_b: usize,
    rules: &KineticRules,
) -> f64 {
    let (ma, mb) = (mass_kg(mass_a_mda), mass_kg(mass_b_mda));
    if ma <= 0.0 || mb <= 0.0 {
        return 0.0;
    }
    let mu = ma * mb / (ma + mb);
    let r = radius_m(atoms_a, rules) + radius_m(atoms_b, rules);
    let sigma = std::f64::consts::PI * r * r;
    rules.steric_factor * sigma * (8.0 * BOLTZMANN / (std::f64::consts::PI * mu)).sqrt()
}

/// The structural part of a **termolecular** prefactor: a bimolecular encounter
/// that must also find a third body nearby.
pub fn termolecular_prefactor(
    mass_a_mda: i64,
    atoms_a: usize,
    mass_b_mda: i64,
    atoms_b: usize,
    rules: &KineticRules,
) -> f64 {
    bimolecular_prefactor(mass_a_mda, atoms_a, mass_b_mda, atoms_b, rules)
        * rules.encounter_volume_m3
}

/// The full Arrhenius prefactor at a temperature, given the stored structural
/// constant and the reaction's molecularity.
///
/// This is the one place that knows how to reconstitute `A(T)`, so the meaning
/// of the stored number and its use can never drift apart.
#[inline]
pub fn prefactor_at(structural: f64, molecularity: u32, temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 {
        return 0.0;
    }
    match molecularity {
        0 | 1 => structural * BOLTZMANN * temperature_k / PLANCK,
        _ => structural * temperature_k.sqrt(),
    }
}
