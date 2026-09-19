//! Self-assembly — when amphiphiles stop being dissolved and start being a
//! surface.
//!
//! # The one thing surfactants do that nothing else does
//!
//! Dissolve salt in water and the dissolved concentration rises until the water
//! is saturated. Dissolve a surfactant and something stranger happens: the free
//! concentration rises to a certain value — the **critical micelle
//! concentration** — and then *stops*, however much more you add. Every further
//! molecule joins an aggregate instead. The free monomer is buffered.
//!
//! That behaviour is the analytic anchor of this module, because it is what
//! makes membranes possible at all. A structure whose concentration self-limits
//! grows by accretion rather than by dilution: add material and you get *more
//! surface*, not a stronger solution. Bags, not brine.
//!
//! # Derived, never stored
//!
//! There is no membrane column and no membrane object. The assembled material in
//! a cell is a **function of the populations already there**, recomputed when
//! needed — exactly as temperature is computed from energy and never stored
//! (Phase 2), and for the same reasons: nothing to keep in sync, nothing to
//! serialize, nothing that can disagree with the state it came from. A world's
//! membranes survive save and resume because they were never separately saved.
//!
//! # Real units, and why that finally works
//!
//! Every number here is physical, and since Phase 5e the chemistry underneath is
//! too, so they connect without anyone tuning anything:
//!
//!   * `cmc_ref_per_m3` defaults to 10²⁴ m⁻³ — about 1.7 mmol/L, where real
//!     short-chain surfactants aggregate.
//!   * `area_per_molecule_m2` defaults to 4×10⁻¹⁹ m² — 0.4 nm², a lipid
//!     headgroup.
//!
//! Those two constants say something concrete about scale. Closing a bag needs
//! `enclosing_area / area_per_molecule` molecules:
//!
//! ```text
//!     1 m³      (a tank)        1.2 × 10¹⁹
//!     10⁻¹⁵ m³  (a eukaryote)   1.2 × 10⁹
//!     10⁻¹⁸ m³  (a bacterium)   1.2 × 10⁷
//! ```
//!
//! That last figure is worth dwelling on, because nothing here was fitted to it:
//! *E. coli*'s membrane holds on the order of 10⁷ lipid molecules. A model given
//! only a headgroup area and the geometry of a sphere lands on the size of a
//! real bacterial membrane. Real cell, real chemistry, real numbers, agreeing —
//! and they only stop agreeing if somebody fudges one of them.
//!
//! The critical concentration itself follows from amphiphilicity
//! ([`crate::amphiphile`]): the stronger the amphiphile, the lower the
//! concentration at which it prefers company to solvent. Real surfactants drop
//! their CMC roughly tenfold for every two carbons added to the tail, and the
//! amphiphilicity score is linear in tail length at a *measured* 3.66 per
//! carbon, so one decade is 7.3 units of score:
//!
//! ```text
//!     cmc(s) = cmc_ref · 10^(−(score − min_score) / cmc_decade)
//! ```
//!
//! `cmc_decade` is therefore not a free parameter but the empirical decade of
//! this model, in this model's units.

use crate::amphiphile::amphiphilicity;
use crate::element::PeriodicTable;
use crate::molecule::Molecule;

/// The declared scales of self-assembly. Every one of them is a physical
/// quantity with a measured value, not a fitting constant.
#[derive(Clone, Copy, Debug)]
pub struct MembraneRules {
    /// Critical concentration, molecules·m⁻³, of an amphiphile whose score is
    /// exactly `min_score`. Default 10²⁴ ≈ 1.7 mmol/L.
    pub cmc_ref_per_m3: f64,
    /// Amphiphilicity score that lowers the critical concentration tenfold.
    /// Traube's rule says two carbons; the score rises 3.66 per carbon, so this
    /// is 7.3.
    pub cmc_decade: f64,
    /// Below this score a molecule never assembles, however concentrated.
    pub min_score: f64,
    /// Area one assembled molecule occupies in the surface, m². Default
    /// 4×10⁻¹⁹ — 0.4 nm², a lipid headgroup.
    pub area_per_molecule_m2: f64,
    /// Permeability of a fully closed compartment. Not zero: real bilayers leak,
    /// which is what lets a protocell take anything up at all.
    pub residual_permeability: f64,
}

impl Default for MembraneRules {
    fn default() -> Self {
        MembraneRules {
            cmc_ref_per_m3: 1.0e24,
            cmc_decade: 7.3,
            min_score: 3.0,
            area_per_molecule_m2: 4.0e-19,
            residual_permeability: 0.01,
        }
    }
}

/// The critical concentration of a species, molecules·m⁻³, or `None` if it is
/// not amphiphilic enough to assemble at any concentration.
pub fn critical_concentration(score: f64, rules: &MembraneRules) -> Option<f64> {
    if !(score >= rules.min_score) {
        return None;
    }
    let decades = (score - rules.min_score) / rules.cmc_decade.max(f64::MIN_POSITIVE);
    Some(rules.cmc_ref_per_m3 * 10f64.powf(-decades))
}

/// How many molecules of one species are surface rather than solution.
///
/// Free monomer is capped at `cmc · volume`; everything beyond that is surface.
/// The split is exact in integers — `free + assembled == pop` always — so this
/// can never invent or lose a molecule, whatever the floats do.
pub fn assembled_count(pop: i128, score: f64, volume_m3: f64, rules: &MembraneRules) -> i128 {
    if pop <= 0 {
        return 0;
    }
    let Some(cmc) = critical_concentration(score, rules) else { return 0 };
    let free_cap = cmc * volume_m3;
    if !(free_cap < i128::MAX as f64) {
        return 0;
    }
    let free = (free_cap.floor() as i128).min(pop).max(0);
    pop - free
}

/// Total assembled material in a cell, in molecules, over every species.
///
/// `scores` is indexed by species id — cached by the caller, since a molecule
/// cannot change once interned and its score is therefore permanent.
pub fn membrane_molecules(
    pops: &[i128],
    scores: &[f64],
    volume_m3: f64,
    rules: &MembraneRules,
) -> i128 {
    let mut total = 0i128;
    for (s, &pop) in pops.iter().enumerate() {
        total += assembled_count(pop, scores.get(s).copied().unwrap_or(0.0), volume_m3, rules);
    }
    total
}

/// The surface area needed to enclose a volume, m² — the area of the sphere of
/// that volume, `(36π)^(1/3) · V^(2/3)`.
///
/// A sphere is the least-area way to enclose a volume, so this is the *minimum*
/// a compartment could get away with, and any real bag needs at least this much.
pub fn enclosing_area(volume_m3: f64) -> f64 {
    if volume_m3 <= 0.0 {
        return 0.0;
    }
    (36.0 * std::f64::consts::PI).cbrt() * volume_m3.powf(2.0 / 3.0)
}

/// What fraction of the boundary the assembled material can cover.
///
/// Below 1 the bag is not closed and things pass through the gaps. At 1 it is
/// just sealed. Above 1 there is more boundary than this volume needs — which is
/// the beginning of a division, since two smaller bags have more total surface
/// than one large one.
pub fn coverage(assembled: i128, volume_m3: f64, rules: &MembraneRules) -> f64 {
    let area = enclosing_area(volume_m3);
    if area <= 0.0 || assembled <= 0 {
        return 0.0;
    }
    (assembled as f64 * rules.area_per_molecule_m2) / area
}

/// How freely a cell exchanges with its neighbours, given its coverage.
///
/// Linear in the uncovered fraction down to a residual: an unbroken bilayer is
/// not a perfect wall, and if it were, nothing could ever get in and a
/// compartment would be a tomb rather than a cell.
pub fn permeability(coverage: f64, rules: &MembraneRules) -> f64 {
    let open = (1.0 - coverage).clamp(0.0, 1.0);
    rules.residual_permeability + (1.0 - rules.residual_permeability) * open
}

/// Amphiphilicity for every molecule in registry order, ready to cache.
pub fn score_table<'a>(
    molecules: impl Iterator<Item = &'a Molecule>,
    table: &PeriodicTable,
) -> Vec<f64> {
    molecules.map(|m| amphiphilicity(m, table)).collect()
}
