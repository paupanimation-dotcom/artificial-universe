//! Catalysis, derived from structure — and the door it opens.
//!
//! # Why this module has to exist
//!
//! Phase 5a's moves — Form, Split, Raise, Lower — join two molecules into one,
//! break one into two, or shuffle a bond order in place. Their stoichiometry is
//! therefore always 2→1, 1→2, or 1→1, and a species never appears on both sides
//! of a derived reaction. (`network_has_no_stoichiometric_autocatalysis` proves
//! this over a real enumerated network, not by assertion.)
//!
//! That is a stronger statement than it first looks. **A molecule that cannot
//! appear among its own reaction's reactants cannot make more of itself.** No
//! amount of running, no temperature, no transport will produce self-replication
//! in such a network, because the shape of the moves forbids it. We went looking
//! for the road to a replicator and found the road did not exist.
//!
//! Catalysis is what builds it. A catalyst enters a reaction and leaves it
//! unchanged — it appears on *both* sides — and the moment the catalyst is also
//! a product, the stoichiometry reads
//!
//! ```text
//!     A + B + AB  →  2 AB
//! ```
//!
//! which is a molecule making more of itself — one copy of AB holding A and B
//! together while a second copy is assembled from them. Autocatalysis is not
//! added here; the *possibility* of it is. Whether any world's chemistry actually realises
//! it is a question for [`crate::raf`] to answer by looking.
//!
//! # The rule, and its honest name: a template model
//!
//! A catalyst works by holding reactants together and stabilising the state
//! between them. This module renders that crudely but structurally:
//!
//! **C catalyses the joining of A and B if C can grip both at once.** Concretely,
//! C must have two free valence *slots* — two spare hands, which may both belong
//! to one bridging atom — one gripping A and one gripping B, where a grip is
//! scored by the same bond-energy model the rest
//! of chemistry uses — so gripping strength is a fact about elements, not a
//! table someone typed. Only *polar* grips count (`min_grip_j_mol`): a site that
//! grips must be chemically unlike the thing it grips, which is the crude form
//! of the true statement that dissimilar atoms attract.
//!
//! What this captures: proximity, orientation, the fact that catalysis is a
//! property of structure and that some molecules are better at it than others.
//! What it does *not* capture: electronic stabilisation of a specific transition
//! state, geometry, strain, solvent, or anything resembling an enzyme's active
//! site. It is a coarse model, stated as one. Its virtue is that no molecule was
//! ever declared to be a catalyst — the relation is computed from graphs the
//! world drew for itself, so a species invented on tick 900 is scored by exactly
//! the same rule as one present at boot.
//!
//! # The law a catalyst may not break
//!
//! **A catalyst changes rates. It never changes equilibria.** This is not a
//! stylistic preference; it is thermodynamics, and it is the analytic anchor
//! this phase is validated against. A catalyst that shifted an equilibrium would
//! be a perpetual-motion machine: run the reaction forward over the catalyst,
//! back without it, and extract work from the cycle forever.
//!
//! The implementation is obliged to obey it, and does so by construction. Phase
//! 5a built barriers as `Ea_f = B + max(ΔH, 0)` and `Ea_r = B + max(−ΔH, 0)`, so
//! `Ea_f − Ea_r = ΔH` exactly and `K = exp(−ΔH/RT)`. Catalysis subtracts the
//! **same absolute amount** from both directions, so the difference — and
//! therefore K — is untouched to the last bit. Two things follow, and both are
//! load-bearing:
//!
//!   * The reduction must be an absolute quantity in J·mol⁻¹, never a fraction
//!     of each barrier. A fixed *fraction* would scale `Ea_f` and `Ea_r` by the
//!     same factor and change their difference — quietly moving every catalysed
//!     equilibrium. This is the kind of error that produces a plausible-looking
//!     simulation which is silently a free-energy machine.
//!   * The reduction is capped strictly below the intrinsic barrier `B`, which
//!     is the smaller of the two barriers. So neither direction can be driven to
//!     zero or negative, and the identity survives every clamp.

use crate::element::{PeriodicTable, Z};
use crate::molecule::Molecule;
use crate::network::NetworkRules;
use crate::reaction::{BondEnergyModel, Reaction, SpeciesId, SpeciesRegistry, Term};

/// The dials of the template model. Like `NetworkRules`, these are the declared
/// *scale* of a derived relation — never a list of which molecule catalyses what.
#[derive(Clone, Copy, Debug)]
pub struct CatalysisRules {
    /// The largest share of the intrinsic barrier `B` a perfect catalyst may
    /// remove. Strictly below 1: the barrier is lowered, never abolished, so
    /// both directions stay positive and `Ea_f − Ea_r = ΔH` survives.
    pub max_reduction_frac: f64,

    /// A grip weaker than this (J·mol⁻¹) does not count at all. Since every bond
    /// is worth at least `base_per_order`, this threshold is what makes the
    /// relation *sparse*: only grips with real polarity qualify. Set it below
    /// `base_per_order` and every molecule with two free sites becomes a
    /// universal catalyst, which is both bad chemistry and a useless signal.
    pub min_grip_j_mol: f64,

    /// The grip strength at which a catalyst achieves half its maximum effect.
    /// Gives a saturating response rather than a cliff.
    pub grip_half_j_mol: f64,
}

impl Default for CatalysisRules {
    fn default() -> Self {
        CatalysisRules {
            max_reduction_frac: 0.75,
            min_grip_j_mol: 250_000.0,
            grip_half_j_mol: 150_000.0,
        }
    }
}

/// Every unused valence **slot** in `mol`, as (atom index, element) — an atom
/// with two spare slots appears twice.
///
/// Counting slots rather than atoms matters, and getting it wrong was a real
/// bug. The first version returned one entry per atom, so a catalyst had to
/// offer its two grips on two *different* atoms. That silently excluded the
/// commonest bridge in chemistry: a single atom holding two partners at once,
/// which is exactly what the oxygen in A–O–B does. In a world of small
/// molecules it excluded nearly every catalyst there was, and the integration
/// tests caught it as "enabling catalysis changed nothing".
fn free_sites(mol: &Molecule, table: &PeriodicTable) -> Vec<(usize, Z)> {
    let mut out = Vec::new();
    for (i, a) in mol.atoms().iter().enumerate() {
        let cap = match table.get(a.z) {
            Some(e) => e.valence,
            None => continue,
        };
        let used = mol.used_valence(i);
        for _ in used..cap {
            out.push((i, a.z));
        }
    }
    out
}

/// The strength of a single transient bond between two elements, scored by the
/// same model that prices every real bond in the engine. Reusing it is the point:
/// a catalyst's grip and a product's stability are quoted in the same currency.
fn grip_strength(model: &BondEnergyModel, table: &PeriodicTable, a: Z, b: Z) -> f64 {
    let (ea, eb) = match (table.get(a), table.get(b)) {
        (Some(x), Some(y)) => (x, y),
        _ => return 0.0,
    };
    let d = (ea.electronegativity() - eb.electronegativity()).abs();
    model.base_per_order + model.electroneg_coefficient * d * d
}

/// The best grip a *specific* catalyst atom can get on any free site of `target`.
fn grip_on(
    model: &BondEnergyModel,
    table: &PeriodicTable,
    catalyst_site: Z,
    target_sites: &[(usize, Z)],
) -> f64 {
    target_sites
        .iter()
        .map(|&(_, z)| grip_strength(model, table, catalyst_site, z))
        .fold(0.0, f64::max)
}

/// How much barrier (J·mol⁻¹) catalyst `c` removes from a reaction that joins or
/// separates `a` and `b` — or `None` if it cannot grip both.
///
/// The catalyst needs two free valence **slots**: one molecule cannot hold two
/// partners with the same hand, but an atom with two spare hands may hold both
/// itself. The hold is the *weaker* of the two
/// grips, because a bridge is only as good as its weaker end, and the response
/// saturates so that an enormously polar site does not buy an unlimited discount.
pub fn catalytic_reduction(
    c: &Molecule,
    a: &Molecule,
    b: &Molecule,
    table: &PeriodicTable,
    model: &BondEnergyModel,
    rules: &CatalysisRules,
    net: &NetworkRules,
) -> Option<f64> {
    let c_sites = free_sites(c, table);
    if c_sites.len() < 2 {
        return None;
    }
    let a_sites = free_sites(a, table);
    let b_sites = free_sites(b, table);
    if a_sites.is_empty() || b_sites.is_empty() {
        return None;
    }

    // Best pair of distinct catalyst atoms, one per partner. Grids this small
    // make the exhaustive search free, and exhaustive means order-independent.
    let mut best_hold = 0.0f64;
    for (pi, &(_, zi)) in c_sites.iter().enumerate() {
        for (pj, &(_, zj)) in c_sites.iter().enumerate() {
            // Two *slots*, which may belong to one bridging atom.
            if pi == pj {
                continue;
            }
            let ga = grip_on(model, table, zi, &a_sites);
            let gb = grip_on(model, table, zj, &b_sites);
            if ga < rules.min_grip_j_mol || gb < rules.min_grip_j_mol {
                continue;
            }
            let hold = ga.min(gb);
            if hold > best_hold {
                best_hold = hold;
            }
        }
    }
    if best_hold <= 0.0 {
        return None;
    }

    // Absolute, capped strictly below the intrinsic barrier: see the module docs
    // on why a *fraction* of each barrier would be a free-energy machine.
    let cap = rules.max_reduction_frac * net.intrinsic_barrier_j_mol;
    let saturating = best_hold / (best_hold + rules.grip_half_j_mol);
    Some(cap * saturating)
}

/// The two molecules a reaction joins or separates, if it is that kind of
/// reaction. Form gives them as reactants, Split as products — either way the
/// pair is the same, which is what makes a catalyst's effect identical in both
/// directions and the equilibrium exactly untouched.
fn joined_pair(r: &Reaction) -> Option<(SpeciesId, SpeciesId)> {
    let side = |terms: &[Term]| -> Option<(SpeciesId, SpeciesId)> {
        match terms {
            [x] if x.count == 2 => Some((x.species, x.species)),
            [x, y] if x.count == 1 && y.count == 1 => Some((x.species, y.species)),
            _ => None,
        }
    };
    if r.reactants.len() == 1 && r.reactants[0].count == 1 {
        side(&r.products)
    } else if r.products.len() == 1 && r.products[0].count == 1 {
        side(&r.reactants)
    } else {
        None
    }
}

/// Every species that catalyses reaction `r`, drawn from `candidates`.
///
/// # Which species are excluded, and why it must be the pair
///
/// A catalyst may not be one of the two fragments being joined or separated —
/// those are reagents, consumed and produced, not catalysts.
///
/// It may, however, be the **joined molecule itself**. That is not a loophole;
/// it is textbook autocatalysis: one copy of AB grips A and B and speeds their
/// joining while a *second* copy is made from them. The catalysing copy leaves
/// unchanged, which is the only thing being a catalyst requires.
///
/// Note carefully that the exclusion is stated in terms of the *pair* (A, B) and
/// not in terms of "appears among this reaction's reactants or products". The
/// first version of this function used the latter, and it was wrong in a way
/// worth recording: for `A + B → AB` the molecule AB is a product, and for the
/// inverse `AB → A + B` it is a reactant, so a sides-based rule would have let
/// AB catalyse the forward reaction and forbidden it from catalysing the
/// reverse. A catalyst present in one direction only is precisely a machine that
/// shifts an equilibrium — the thermodynamic sin this module's own documentation
/// forbids. Excluding the pair is symmetric under inversion by construction,
/// because [`joined_pair`] returns the same pair for both directions.
pub fn catalysts_for(
    r: &Reaction,
    candidates: &[SpeciesId],
    reg: &SpeciesRegistry,
    table: &PeriodicTable,
    model: &BondEnergyModel,
    rules: &CatalysisRules,
    net: &NetworkRules,
) -> Vec<(SpeciesId, f64)> {
    let (sa, sb) = match joined_pair(r) {
        Some(p) => p,
        None => return Vec::new(),
    };
    let (a, b) = match (reg.get(sa), reg.get(sb)) {
        (Some(x), Some(y)) => (x, y),
        _ => return Vec::new(),
    };
    let mut out = Vec::new();
    for &cid in candidates {
        if cid == sa || cid == sb {
            continue;
        }
        let c = match reg.get(cid) {
            Some(m) => m,
            None => continue,
        };
        if let Some(red) = catalytic_reduction(c, a, b, table, model, rules, net) {
            out.push((cid, red));
        }
    }
    out.sort_by(|x, y| x.0 .0.cmp(&y.0 .0));
    out
}

/// Merge duplicate species into a single term, ascending by id.
///
/// Needed because a catalyst can already appear on a side: catalysing the split
/// `AB → A + B` with AB gives reactants `[AB, AB]`, which is exactly `2 AB` and
/// should be written that way. Sorted output keeps the reaction's canonical key
/// stable, so the dedup in `enumerate_reactions` still sees one reaction.
fn coalesce(mut terms: Vec<Term>) -> Vec<Term> {
    terms.sort_by_key(|t| t.species.0);
    let mut out: Vec<Term> = Vec::with_capacity(terms.len());
    for t in terms {
        match out.last_mut() {
            Some(last) if last.species == t.species => last.count += t.count,
            _ => out.push(t),
        }
    }
    out
}

/// A catalysed reaction, ready for the reactor.
#[derive(Clone, Debug)]
pub struct Catalysed {
    /// The reaction as the reactor should run it: catalyst on both sides,
    /// barrier lowered in both directions.
    pub reaction: Reaction,
    /// Which species is doing the catalysing.
    pub catalyst: SpeciesId,
    /// Index of the uncatalysed parent in the input slice.
    pub parent: usize,
}

/// Derive the catalysed variants of a reaction network.
///
/// The uncatalysed reactions are *kept* by the caller: a catalyst makes a road
/// faster, it does not close the old one. What is returned here is the extra
/// roads, each one a mass-action reaction with the catalyst as both reactant and
/// product — so its rate is proportional to the catalyst's concentration, which
/// is what makes catalysis a dynamical thing rather than a bookkeeping trick.
///
/// Note the shape of the output when a catalyst is also a product of its own
/// parent reaction: `A + B + AB → 2 AB`. Nothing here special-cases that. It is
/// simply what the general rule says when the structure happens to line up, and
/// it is the first time in this engine's life that a species can increase its
/// own number.
pub fn catalysed_variants(
    reactions: &[Reaction],
    candidates: &[SpeciesId],
    reg: &SpeciesRegistry,
    table: &PeriodicTable,
    model: &BondEnergyModel,
    rules: &CatalysisRules,
    net: &NetworkRules,
) -> Vec<Catalysed> {
    let mut out = Vec::new();
    for (i, r) in reactions.iter().enumerate() {
        for (cid, reduction) in catalysts_for(r, candidates, reg, table, model, rules, net) {
            let mut cr = r.clone();
            // The same absolute reduction in both directions — the forward and
            // reverse variants are generated from the same pair, so `Ea_f − Ea_r`
            // and therefore K are preserved exactly.
            cr.activation_j = (cr.activation_j - reduction).max(0.0);
            cr.reactants.push(Term { species: cid, count: 1 });
            cr.products.push(Term { species: cid, count: 1 });
            cr.reactants = coalesce(std::mem::take(&mut cr.reactants));
            cr.products = coalesce(std::mem::take(&mut cr.products));
            // The variant's prefactor is *recomputed* from its own reactant
            // set, never scaled from the parent's. Molecularity has gone up by
            // one, and a prefactor's units depend on molecularity — so scaling
            // a unimolecular parent's dimensionless steric factor by an
            // encounter volume would not yield a bimolecular prefactor at all.
            // (It did, briefly, and produced a crossover concentration ten
            // orders of magnitude denser than solid matter, which is how the
            // test caught it.)
            //
            // Either way the catalyst must be *found* as well as the substrates,
            // and that is the physical reason behind Phase 5c's empirical
            // discovery that catalysis only pays when the barrier removed is
            // large.
            cr.pre_exponential =
                crate::network::structural_prefactor(&cr.reactants, reg, table, &net.kinetics);
            out.push(Catalysed { reaction: cr, catalyst: cid, parent: i });
        }
    }
    out
}
