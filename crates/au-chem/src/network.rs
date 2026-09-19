//! The generative reaction network — chemistry that discovers.
//!
//! # Why this module exists
//!
//! Everything before this point ran on a *declared* chemistry: config listed the
//! molecules that could exist and the reactions that could fire, the way it lists
//! materials. That was the honest Phase 3 scope, and it was also a ceiling — a
//! world whose molecule list is written down in advance can never contain a
//! surprise, and the whole point of the layers above chemistry is that they must
//! be able to. **Life cannot emerge from a vocabulary that is closed.**
//!
//! This module opens it. Reactions are no longer declared; they are *derived from
//! molecular structure* by a small set of local graph moves, and when a move
//! produces a molecule the world has never seen, that molecule is interned on the
//! spot — a new species, discovered, not designed. Put hydrogen and oxygen atoms
//! in a hot cell and nobody has to tell the engine about water: bonds form
//! because forming them is possible, and water appears because it is reachable
//! and — being the deepest energy well in the neighbourhood — it stays.
//!
//! # The moves, and why they come in mutually inverse pairs
//!
//! A candidate reaction is one application of one move:
//!
//!   * **Form** — a new single bond between a free-valence atom of one molecule
//!     and a free-valence atom of another: `A + B → A–B`.
//!   * **Split** — remove a single bond whose loss disconnects the molecule:
//!     `A–B → A + B`.
//!   * **Raise** — increase an existing bond's order by one, where both atoms
//!     have a free valence slot to pay with: `A–B → A=B`, `A=B → A≡B`.
//!   * **Lower** — the opposite: `A≡B → A=B`, `A=B → A–B`.
//!
//! Form/Split are exact inverses, and so are Raise/Lower. That pairing is not
//! tidiness — it is thermodynamics. Every reaction the generator can produce, it
//! can also produce the reverse of (once the products exist to react), so every
//! transformation is reversible and the network can settle into genuine
//! equilibria instead of one-way ratchets. The activation model below then
//! guarantees those equilibria sit exactly where Boltzmann says they should.
//!
//! Deliberately **not** generated, each a documented exclusion rather than an
//! oversight:
//!
//!   * *Ring closure* (a bond between two atoms of the same molecule) and its
//!     inverse, ring opening (a split that leaves the molecule connected). Rings
//!     matter enormously later — aromatics, sugars — but admitting them means
//!     admitting intramolecular moves as a pair, and Phase 5a's validation world
//!     (H and O) does not need them. When they come, they come together.
//!   * *Concerted substitutions* (A–B + C → A–C + B in one step). Real, but
//!     representable as split-then-form through an intermediate; collapsing them
//!     into one move changes kinetics, not reachability.
//!
//! # The activation model, and the one identity it must satisfy
//!
//! A generated reaction needs an activation energy, and the model is deliberately
//! the simplest one that preserves the thermodynamics the engine already
//! validated:
//!
//! ```text
//!     Ea_forward = B + max(ΔH, 0)
//!     Ea_reverse = B + max(−ΔH, 0)
//! ```
//!
//! where `B` is a single declared intrinsic barrier and ΔH is the reaction
//! enthalpy the bond model derives. Subtracting the two lines gives the identity
//! everything rests on:
//!
//! ```text
//!     Ea_forward − Ea_reverse = ΔH        (exactly, by construction)
//! ```
//!
//! With equal pre-exponentials both ways, the equilibrium constant is then
//! `K = exp(−ΔH/RT)` — Boltzmann — and van 't Hoff behaviour follows for free.
//! Phase 3 *tested* that identity on hand-declared reactions; here it is
//! *structural*, true of every reaction the generator can ever emit. A test
//! still pins it, because "true by construction" is a claim about code that
//! refactoring can silently break.
//!
//! The per-event enthalpy stored on the reaction is the molar enthalpy divided by
//! Avogadro's number — the honest bridge from the bond model's J/mol to a cell
//! that counts individual molecules. (Activation energies stay molar, because the
//! Arrhenius exponent divides by R·T.)
//!
//! # Symmetry, and why enumeration does not explode
//!
//! Attaching a hydrogen to any of benzene's six carbons gives the same product,
//! and enumerating it six times would multiply candidates for nothing. The
//! generator therefore works on **symmetry classes** (`Molecule::symmetry_classes`
//! — the converged Morgan ranks): one representative site per class for Form, one
//! representative bond per class for Split/Raise/Lower. Water offers one O–H bond
//! to split, not two.
//!
//! Refinement can, on certain regular graphs, merge classes it should not — so
//! finished candidates are *also* deduplicated by the canonical forms of their
//! full reactant and product lists. The class machinery keeps enumeration small;
//! the canonical dedup keeps it correct. Neither is trusted alone.
//!
//! # The honest cap
//!
//! `max_atoms` bounds the size of any molecule the generator will create. This is
//! the combinatorial cliff guard, the chemistry analogue of the chunk budget: an
//! unbounded generator on an exothermic polymerising chemistry will happily build
//! molecules the size of the heap. Reactions whose product would exceed the cap
//! are simply not offered — a declared edge of the world, stated rather than hit.

use crate::element::PeriodicTable;
use crate::molecule::{Bond, BondOrder, Molecule};
use crate::kinetics::{
    bimolecular_prefactor, termolecular_prefactor, unimolecular_prefactor, KineticRules,
};
use crate::reaction::{BondEnergyModel, Reaction, SpeciesId, SpeciesRegistry, Term};
use au_physics::Energy;
use std::collections::{BTreeMap, BTreeSet};

/// Avogadro's number — the bridge between the bond model's molar energies and a
/// cell that counts molecules one by one.
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// The declared parameters of an open chemistry. Everything else is derived.
#[derive(Clone, Copy, Debug)]
pub struct NetworkRules {
    /// The intrinsic activation barrier `B`, J/mol — the cost of *any*
    /// rearrangement, over and above the thermodynamic hill. One number for the
    /// whole chemistry; making it structural (per bond type) is a later
    /// refinement that changes kinetics, never equilibria.
    pub intrinsic_barrier_j_mol: f64,
    /// The declared scales of collision kinetics. The prefactor itself is no
    /// longer declared: it is computed per reaction from the reactants' masses
    /// and sizes (see [`crate::kinetics`]).
    ///
    /// Note what this changes about equilibrium. Forward and reverse prefactors
    /// used to be equal by construction, which made `K = exp(−ΔH/RT)` exact.
    /// They are no longer equal — an association is bimolecular forwards and
    /// unimolecular backwards — and their ratio is precisely the reaction
    /// entropy, so `K = (A_f/A_r)·exp(−ΔH/RT)`. The ΔS term that was missing for
    /// three phases is not added anywhere; it falls out of counting molecules on
    /// each side of the arrow.
    pub kinetics: KineticRules,
    /// No generated molecule may exceed this many atoms. The cliff guard.
    pub max_atoms: usize,
}

impl Default for NetworkRules {
    fn default() -> Self {
        NetworkRules {
            intrinsic_barrier_j_mol: 60_000.0,
            kinetics: KineticRules::default(),
            max_atoms: 8,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Graph moves — pure functions from molecules to molecules.
// ─────────────────────────────────────────────────────────────────────────────

/// Join two molecules with a new single bond between atom `ai` of `a` and atom
/// `bi` of `b`. The caller has checked valence; this only does the bookkeeping.
pub fn join(a: &Molecule, ai: usize, b: &Molecule, bi: usize) -> Molecule {
    let mut atoms = a.atoms().to_vec();
    atoms.extend_from_slice(b.atoms());
    let off = a.atom_count() as u16;
    let mut bonds = a.bonds().to_vec();
    bonds.extend(b.bonds().iter().map(|d| Bond::new(d.a + off, d.b + off, d.order)));
    bonds.push(Bond::new(ai as u16, bi as u16 + off, BondOrder::Single));
    Molecule::new(atoms, bonds)
}

/// Remove bond `bond_idx` from `m`. If the molecule falls into two pieces,
/// return them; if it stays connected (the bond was part of a ring), return
/// `None` — ring opening is an excluded move, see the module docs.
pub fn split(m: &Molecule, bond_idx: usize) -> Option<(Molecule, Molecule)> {
    let n = m.atom_count();
    // Flood-fill from atom 0 over every bond except the removed one.
    let mut comp = vec![u8::MAX; n];
    let mut stack = vec![0usize];
    comp[0] = 0;
    while let Some(i) = stack.pop() {
        for (k, b) in m.bonds().iter().enumerate() {
            if k == bond_idx {
                continue;
            }
            let j = if b.a as usize == i {
                b.b as usize
            } else if b.b as usize == i {
                b.a as usize
            } else {
                continue;
            };
            if comp[j] == u8::MAX {
                comp[j] = 0;
                stack.push(j);
            }
        }
    }
    if comp.iter().all(|&c| c == 0) {
        return None; // still connected: a ring bond
    }
    // Everything not reached is the second fragment.
    for c in comp.iter_mut() {
        if *c == u8::MAX {
            *c = 1;
        }
    }
    let build = |which: u8| -> Molecule {
        let mut map = vec![u16::MAX; n];
        let mut atoms = Vec::new();
        for i in 0..n {
            if comp[i] == which {
                map[i] = atoms.len() as u16;
                atoms.push(m.atoms()[i]);
            }
        }
        let bonds = m
            .bonds()
            .iter()
            .enumerate()
            .filter(|(k, b)| *k != bond_idx && comp[b.a as usize] == which)
            .map(|(_, b)| Bond::new(map[b.a as usize], map[b.b as usize], b.order))
            .collect();
        Molecule::new(atoms, bonds)
    };
    Some((build(0), build(1)))
}

/// The same molecule with bond `bond_idx` at a different order.
pub fn with_order(m: &Molecule, bond_idx: usize, order: BondOrder) -> Molecule {
    let mut bonds = m.bonds().to_vec();
    bonds[bond_idx] = Bond::new(bonds[bond_idx].a, bonds[bond_idx].b, order);
    Molecule::new(m.atoms().to_vec(), bonds)
}

// ─────────────────────────────────────────────────────────────────────────────
// Site enumeration — one representative per symmetry class.
// ─────────────────────────────────────────────────────────────────────────────

/// One representative atom index per symmetry class that still has free valence.
fn free_sites(m: &Molecule, table: &PeriodicTable) -> Vec<usize> {
    let classes = m.symmetry_classes();
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for i in 0..m.atom_count() {
        let Some(el) = table.get(m.atoms()[i].z) else { continue };
        if m.used_valence(i) >= el.valence {
            continue;
        }
        if seen.insert(classes[i]) {
            out.push(i);
        }
    }
    out
}

/// One representative bond index per bond symmetry class. A bond's class is the
/// unordered pair of its endpoints' classes plus its order — two bonds in the
/// same class connect indistinguishable environments, so acting on either yields
/// isomorphic products.
fn bond_classes(m: &Molecule) -> Vec<usize> {
    let classes = m.symmetry_classes();
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for (k, b) in m.bonds().iter().enumerate() {
        let (ca, cb) = (classes[b.a as usize], classes[b.b as usize]);
        let key = (ca.min(cb), ca.max(cb), b.order.slots());
        if seen.insert(key) {
            out.push(k);
        }
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// The generator.
// ─────────────────────────────────────────────────────────────────────────────

/// A normalized multiset of terms, used as half a dedup key.
fn term_key(terms: &[Term], reg: &SpeciesRegistry) -> Vec<u8> {
    let mut parts: Vec<Vec<u8>> = terms
        .iter()
        .map(|t| {
            let mut v = reg.canonical(t.species).unwrap().0.clone();
            v.extend_from_slice(&t.count.to_le_bytes());
            v
        })
        .collect();
    parts.sort();
    let mut out = Vec::new();
    for p in parts {
        out.extend_from_slice(&(p.len() as u32).to_le_bytes());
        out.extend_from_slice(&p);
    }
    out
}

/// Build a reaction from reactant/product terms, deriving its energetics.
fn make_reaction(
    reactants: Vec<Term>,
    products: Vec<Term>,
    reg: &SpeciesRegistry,
    table: &PeriodicTable,
    model: &BondEnergyModel,
    rules: &NetworkRules,
) -> Reaction {
    let reactants_for_kinetics = reactants.clone();
    let dh_molar = model.reaction_enthalpy_molar(&reactants, &products, reg, table);
    let ea = rules.intrinsic_barrier_j_mol + dh_molar.max(0.0);
    let dh_event = dh_molar / AVOGADRO;
    Reaction {
        reactants,
        products,
        activation_j: ea,
        enthalpy: Energy::from_joules(dh_event),
        enthalpy_j: dh_event,
        pre_exponential: structural_prefactor(&reactants_for_kinetics, reg, table, &rules.kinetics),
    }
}

/// The structural prefactor for a reaction, from the molecules that must meet.
///
/// Unimolecular reactions get a steric factor and pick up `k_B·T/h` at rate
/// time; bimolecular ones get a collision cross-section and reduced mass read
/// off the two graphs; anything higher pays an encounter volume on top. All of
/// it comes from the molecules themselves, so a species the world invented
/// yesterday is priced by the same rule as one present at boot.
pub fn structural_prefactor(
    reactants: &[Term],
    reg: &SpeciesRegistry,
    table: &PeriodicTable,
    kinetics: &KineticRules,
) -> f64 {
    let molecularity: u32 = reactants.iter().map(|t| t.count).sum();
    if molecularity <= 1 {
        return unimolecular_prefactor(kinetics);
    }
    // The two partners that actually collide. A term of count two is one
    // species meeting itself.
    let mut partners: Vec<SpeciesId> = Vec::with_capacity(2);
    for t in reactants {
        for _ in 0..t.count {
            if partners.len() < 2 {
                partners.push(t.species);
            }
        }
    }
    if partners.len() < 2 {
        return unimolecular_prefactor(kinetics);
    }
    let describe = |s: SpeciesId| -> (i64, usize) {
        match reg.get(s) {
            Some(m) => (m.mass_mda(table) as i64, m.atom_count()),
            None => (1, 1),
        }
    };
    let (ma, na) = describe(partners[0]);
    let (mb, nb) = describe(partners[1]);
    if molecularity == 2 {
        bimolecular_prefactor(ma, na, mb, nb, kinetics)
    } else {
        termolecular_prefactor(ma, na, mb, nb, kinetics)
    }
}

/// Enumerate every reaction one move away from the given species, interning any
/// product molecule the registry has never seen. **This is where discovery
/// happens** — the registry after this call may be larger than before it.
///
/// Determinism: `present` is sorted and deduplicated internally; sites and bonds
/// are visited in ascending index with class dedup; finished reactions are keyed
/// and emitted in `BTreeMap` order. Two calls from identical state produce
/// identical reactions in identical order, which is what lets discovery live
/// inside a world that must hash the same on every run.
pub fn enumerate_reactions(
    present: &[SpeciesId],
    reg: &mut SpeciesRegistry,
    table: &PeriodicTable,
    model: &BondEnergyModel,
    rules: &NetworkRules,
) -> Vec<Reaction> {
    let mut present: Vec<SpeciesId> = present.to_vec();
    present.sort_by_key(|s| s.0);
    present.dedup();

    let mut out: BTreeMap<Vec<u8>, Reaction> = BTreeMap::new();
    let add = |reactants: Vec<Term>,
                   products: Vec<Term>,
                   reg: &SpeciesRegistry,
                   out: &mut BTreeMap<Vec<u8>, Reaction>| {
        let mut key = term_key(&reactants, reg);
        key.push(0xFF);
        key.extend(term_key(&products, reg));
        out.entry(key)
            .or_insert_with(|| make_reaction(reactants, products, reg, table, model, rules));
    };

    // ── Form: A + B → A–B, one candidate per (site class of A) × (site class of B).
    for ii in 0..present.len() {
        for jj in ii..present.len() {
            let (sa, sb) = (present[ii], present[jj]);
            let (a, b) = (reg.get(sa).unwrap().clone(), reg.get(sb).unwrap().clone());
            if a.atom_count() + b.atom_count() > rules.max_atoms {
                continue;
            }
            for &ai in &free_sites(&a, table) {
                for &bi in &free_sites(&b, table) {
                    // Same species joining itself: sites are interchangeable
                    // across the two copies, so only take ai ≤ bi.
                    if sa == sb && bi < ai {
                        continue;
                    }
                    let product = join(&a, ai, &b, bi);
                    let pid = reg.intern(product);
                    let reactants = if sa == sb {
                        vec![Term { species: sa, count: 2 }]
                    } else {
                        vec![Term { species: sa, count: 1 }, Term { species: sb, count: 1 }]
                    };
                    add(reactants, vec![Term { species: pid, count: 1 }], reg, &mut out);
                }
            }
        }
    }

    // ── Split, Raise, Lower: one candidate per bond class of each present species.
    for &sid in &present {
        let m = reg.get(sid).unwrap().clone();
        for k in bond_classes(&m) {
            let bond = m.bonds()[k];

            // Split — single bonds only (a double bond "splitting" would strand
            // valence; lowering it first is the generated path down).
            if bond.order.slots() == 1 {
                if let Some((f0, f1)) = split(&m, k) {
                    let i0 = reg.intern(f0);
                    let i1 = reg.intern(f1);
                    let products = if i0 == i1 {
                        vec![Term { species: i0, count: 2 }]
                    } else {
                        vec![Term { species: i0, count: 1 }, Term { species: i1, count: 1 }]
                    };
                    add(vec![Term { species: sid, count: 1 }], products, reg, &mut out);
                }
            }

            // Raise — both endpoints must have a free valence slot to pay with.
            let next = BondOrder::from_u8(bond.order.slots() + 1);
            if let Some(next) = next {
                let (va, vb) = (
                    table.get(m.atoms()[bond.a as usize].z).map(|e| e.valence).unwrap_or(0),
                    table.get(m.atoms()[bond.b as usize].z).map(|e| e.valence).unwrap_or(0),
                );
                if m.used_valence(bond.a as usize) < va && m.used_valence(bond.b as usize) < vb {
                    let pid = reg.intern(with_order(&m, k, next));
                    add(
                        vec![Term { species: sid, count: 1 }],
                        vec![Term { species: pid, count: 1 }],
                        reg,
                        &mut out,
                    );
                }
            }

            // Lower — always possible above Single; Single's "lower" is Split.
            if bond.order.slots() > 1 {
                let lower = BondOrder::from_u8(bond.order.slots() - 1).unwrap();
                let pid = reg.intern(with_order(&m, k, lower));
                add(
                    vec![Term { species: sid, count: 1 }],
                    vec![Term { species: pid, count: 1 }],
                    reg,
                    &mut out,
                );
            }
        }
    }

    out.into_values().collect()
}

/// A human-readable formula, "H2O" style, element order as declared in the
/// table. For demos and error messages; nothing in the engine reads it back.
pub fn formula_string(m: &Molecule, table: &PeriodicTable) -> String {
    let counts = m.formula(table);
    let mut s = String::new();
    for (i, el) in table.iter().enumerate() {
        let c = counts[i];
        if c > 0 {
            s.push_str(el.symbol);
            if c > 1 {
                s.push_str(&c.to_string());
            }
        }
    }
    s
}
