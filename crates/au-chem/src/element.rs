//! Elements, as exactly conserved integer counts.
//!
//! # The one decision this whole phase rests on
//!
//! An element is not a substance with properties chosen by a designer. It is a
//! *count of atoms of a given atomic number*, and atoms are neither created nor
//! destroyed by chemistry — only rearranged. Nuclear reactions, which is the only
//! thing that would change an atom's identity, are out of scope. So element
//! counts are conserved **exactly and forever**, by construction, in precisely
//! the way energy and mass already are.
//!
//! This matters for the same reason exact energy mattered in Phase 2, and the
//! reason is worth restating because it is the reason the entire tower stands:
//! **evolution is an adversarial optimiser, and its search space includes the
//! bugs in your chemistry.** If a reaction can quietly produce a carbon atom from
//! nothing, then a self-replicating molecule that exploits that leak will
//! eventually appear — not because anything is clever, but because a free source
//! of a limiting element is the most valuable thing in any chemical world, and
//! selection is an exhaustive search. Conservation cannot be "accurate". It has
//! to be a bit-exact identity, so that there is *no* leak to find.
//!
//! # Why the periodic table is small, and open
//!
//! Real chemistry is built from ~92 natural elements, but life on Earth is built
//! from about six (C, H, O, N, P, S), and the chemistry that *matters* for
//! emergence is overwhelmingly the chemistry of those six plus a handful of ions.
//! So this file does not hardcode a periodic table — it provides the *type* of an
//! element (its conserved count, its mass, its bonding capacity) and lets a world
//! declare which elements exist. A world with an exotic element set is a different
//! universe, which is exactly what PROJECT_VISION asks for.
//!
//! The properties here — valence, electronegativity — are the minimum needed to
//! decide *whether atoms bond and how many bonds they form*. They are deliberately
//! a caricature of real quantum chemistry, because real quantum chemistry is not
//! computable at a million cells (see the crate docs). The bet is that this
//! caricature preserves enough of the causal structure for interesting chemistry
//! to emerge. Whether it does is what the tests, and ultimately Phase 5, decide.

use au_core::hash::{Hasher, WorldHash};

/// Atomic number. The identity of an element — the count of protons — and the one
/// thing chemistry can never change. Two atoms with the same `Z` are the same
/// element; nothing a reaction does alters it.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub struct Z(pub u8);

impl Z {
    pub const H: Z = Z(1);
    pub const C: Z = Z(6);
    pub const N: Z = Z(7);
    pub const O: Z = Z(8);
    pub const P: Z = Z(15);
    pub const S: Z = Z(16);
}

/// The fixed properties of an element.
///
/// "Fixed" is the point: unlike a Phase 2 `Material` (whose bulk properties will,
/// in the fullness of the project, be *derived* from which molecules are present),
/// these are the genuinely atomic facts. An element's mass and proton count do
/// not emerge from anything below — they are the floor.
#[derive(Clone, Debug)]
pub struct Element {
    pub z: Z,
    pub symbol: &'static str,

    /// Atomic mass, in unified atomic mass units (daltons) × 1000, as an integer.
    ///
    /// Integer for the same reason everything conserved is an integer: so that
    /// the total mass of a reaction's inputs equals the total mass of its outputs
    /// *exactly*, with no floating-point residue for a hungry lineage to farm.
    /// ×1000 resolves the isotopic averaging real atomic weights carry (H is
    /// 1.008, not 1).
    pub mass_mda: u32,

    /// Typical number of covalent bonds this atom forms — its valence.
    ///
    /// This is the single most load-bearing simplification in the crate. Real
    /// bonding is a continuous quantum-mechanical affair; here it is an integer
    /// budget. Carbon's valence of 4 is *why* carbon chains, rings, and the entire
    /// combinatorial explosion of organic chemistry are possible — and it is why
    /// this abstraction has any hope of producing life-like complexity. An atom
    /// with all its bonds used is saturated and inert; an atom with a free bond is
    /// reactive. That single rule drives most of what follows.
    pub valence: u8,

    /// Pauling electronegativity × 100, as an integer. How greedily this atom pulls
    /// shared electrons. The *difference* in electronegativity across a bond
    /// decides the bond's polarity and much of its energy — it is why O–H bonds are
    /// strong and polar (and water is water) while C–C bonds are not. Used to
    /// derive bond energies rather than tabulating every pair by hand.
    pub electronegativity_c: u16,
}

impl Element {
    pub fn mass_da(&self) -> f64 {
        self.mass_mda as f64 / 1000.0
    }
    pub fn electronegativity(&self) -> f64 {
        self.electronegativity_c as f64 / 100.0
    }
}

impl WorldHash for Element {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u8(self.z.0);
        h.write_u32(self.mass_mda);
        h.write_u8(self.valence);
        h.write_u16(self.electronegativity_c);
    }
}

/// The set of elements that exist in a given universe.
///
/// Empty by default. A world declares its own chemistry — there is no built-in
/// periodic table, because "which elements exist" is a property of the universe,
/// not of the engine. A validation fixture (`data/reference_chemistry.kv`, loaded
/// only by tests) supplies real measured values so the engine can be checked
/// against known chemistry; it is a measuring instrument, not content.
#[derive(Clone, Debug, Default)]
pub struct PeriodicTable {
    elements: Vec<Element>,
}

impl PeriodicTable {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, e: Element) {
        assert!(
            !self.elements.iter().any(|x| x.z == e.z),
            "element Z={} already declared",
            e.z.0
        );
        self.elements.push(e);
        // Keep sorted by Z so iteration order is deterministic and independent of
        // declaration order — the world hash must not depend on the order a config
        // file happened to list its elements in.
        self.elements.sort_by_key(|e| e.z.0);
    }

    pub fn get(&self, z: Z) -> Option<&Element> {
        self.elements.iter().find(|e| e.z == z)
    }

    pub fn by_symbol(&self, sym: &str) -> Option<&Element> {
        self.elements.iter().find(|e| e.symbol == sym)
    }

    /// The dense index of an element, for compact per-cell storage: a world with
    /// six elements stores six integers per cell, not a sparse map keyed by Z.
    pub fn index_of(&self, z: Z) -> Option<usize> {
        self.elements.iter().position(|e| e.z == z)
    }

    pub fn z_at(&self, index: usize) -> Option<Z> {
        self.elements.get(index).map(|e| e.z)
    }

    pub fn len(&self) -> usize {
        self.elements.len()
    }
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }
    pub fn iter(&self) -> impl Iterator<Item = &Element> {
        self.elements.iter()
    }
}

impl WorldHash for PeriodicTable {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u64(self.elements.len() as u64);
        for e in &self.elements {
            e.hash_into(h);
        }
    }
}

/// A conserved inventory of atoms — a count per element, in one place (a cell, or
/// a molecule population).
///
/// This is the chemical analogue of the Phase 2 energy ledger, and it plays the
/// same role: it is the thing that must balance. Every reaction is a permutation
/// of atoms within this inventory. Nothing enters, nothing leaves, unless it
/// crosses a boundary and is booked — exactly as with energy.
///
/// Counts are `i128`, absurdly large, because a single cell of gas can hold ~10²³
/// molecules and we would rather never think about overflow again. A count is
/// literally a number of atoms; there is no unit and no rounding.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct AtomCount {
    /// Indexed to match a `PeriodicTable` — `counts[i]` is the number of atoms of
    /// the element at index `i`. Dense and cache-friendly.
    counts: Vec<i128>,
}

impl AtomCount {
    pub fn zeros(n_elements: usize) -> Self {
        AtomCount { counts: vec![0; n_elements] }
    }

    pub fn from_counts(counts: Vec<i128>) -> Self {
        AtomCount { counts }
    }

    #[inline]
    pub fn get(&self, index: usize) -> i128 {
        self.counts.get(index).copied().unwrap_or(0)
    }

    #[inline]
    pub fn set(&mut self, index: usize, n: i128) {
        self.counts[index] = n;
    }

    #[inline]
    pub fn add(&mut self, index: usize, n: i128) {
        self.counts[index] = self.counts[index].checked_add(n).expect("atom count overflow");
    }

    /// Move `n` atoms of one element from `self` to `other`. The atomic operation
    /// of chemistry: a symmetric transfer that cannot create or destroy an atom,
    /// only relocate it. Conservation is an identity, not a hope — exactly the
    /// flux-form discipline that got Phase 2 energy to zero error.
    #[inline]
    pub fn transfer_to(&mut self, other: &mut AtomCount, index: usize, n: i128) {
        self.counts[index] = self.counts[index].checked_sub(n).expect("atom underflow");
        other.counts[index] = other.counts[index].checked_add(n).expect("atom overflow");
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }
    pub fn is_empty(&self) -> bool {
        self.counts.iter().all(|&c| c == 0)
    }

    /// Total atoms of every kind. The single number a conservation test watches.
    pub fn total(&self) -> i128 {
        self.counts.iter().sum()
    }

    pub fn as_slice(&self) -> &[i128] {
        &self.counts
    }

    /// Total mass of this inventory, in milli-daltons. Exact: Σ count·mass. Because
    /// both factors are integers, the mass of a reaction balances to the
    /// milli-dalton, which is how "mass is conserved" stops being a slogan and
    /// becomes a checked invariant.
    pub fn mass_mda(&self, table: &PeriodicTable) -> i128 {
        self.counts
            .iter()
            .enumerate()
            .map(|(i, &c)| c * table.iter().nth(i).map(|e| e.mass_mda as i128).unwrap_or(0))
            .sum()
    }
}

impl WorldHash for AtomCount {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u64(self.counts.len() as u64);
        for &c in &self.counts {
            h.write_bytes(&c.to_le_bytes());
        }
    }
}
