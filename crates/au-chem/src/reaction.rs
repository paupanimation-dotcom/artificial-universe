//! Reactions — where chemistry and thermodynamics become one loop.
//!
//! # What a reaction is here
//!
//! A reaction is a rule: *these reactant species, in these amounts, become these
//! product species*. It fires at a rate that depends on how much reactant is
//! present (mass-action kinetics) and how hot it is (Arrhenius temperature
//! dependence) — the same two dependencies that govern every real reaction. And it
//! carries an **energy of reaction** that flows straight into the Phase 2 energy
//! field.
//!
//! That last clause is the whole point of building chemistry *on* physics rather
//! than beside it. An exothermic reaction releases energy into the cell it occurs
//! in, warming it; the warming raises the Arrhenius rate of every reaction there;
//! and an endothermic reaction cools its cell, slowing things down. Heat drives
//! chemistry and chemistry drives heat, through a single shared, conserved energy
//! account. There is no separate "chemical energy" bolted on — there is *energy*,
//! and reactions move it between the bond-configuration ledger and the thermal
//! ledger.
//!
//! # Conservation, twice over, both exact
//!
//! Every reaction must balance atoms exactly (a rule that produces a carbon from
//! nothing is rejected at construction, not at runtime) and must conserve total
//! energy exactly: the energy released as heat is precisely the difference in
//! stored bond energy between reactants and products. Both are integer identities,
//! for the reason that has driven every decision in this project — evolution is an
//! adversarial optimiser, and a chemical world offers it two prizes to farm, free
//! atoms and free energy. Neither may exist.
//!
//! # Why rates, and not just equilibria
//!
//! It would be simpler to jump each reaction straight to its equilibrium. But
//! life is not an equilibrium — it is a system held *far from* equilibrium by a
//! continuous flow of energy and matter (a dissipative structure, exactly like the
//! convection rolls of Phase 2b, exactly like a cell). Kinetics — the fact that
//! reactions take *time*, and that a fast reaction can outrun a
//! thermodynamically-favoured slow one — is what lets metabolism exist: a cell
//! stays alive precisely by kinetically routing matter through pathways that would
//! never dominate at equilibrium. So reactions here have rates, and equilibrium is
//! something the system *approaches* if left alone, never something imposed.

use au_physics::Energy;

use crate::element::PeriodicTable;
use crate::molecule::{CanonicalForm, Molecule};

/// The universal gas constant, in J·mol⁻¹·K⁻¹ × 1000 for a little integer-friendly
/// headroom where it helps. Used in the Arrhenius and van 't Hoff relations.
pub const R_GAS: f64 = 8.314_462_618;

/// The registry of molecular species a world has encountered.
///
/// Species are **discovered, then remembered** — this is the concrete form of "no
/// hardcoded molecules". When a reaction produces a graph, its canonical form is
/// looked up here; if it is new, it is assigned the next `SpeciesId` and recorded.
/// A cell then stores a compact count per `SpeciesId` rather than a pile of graphs.
///
/// The registry only ever grows, and it grows deterministically (species are
/// numbered in discovery order, and discovery order is fixed by the simulation's
/// determinism). Two identical runs assign identical ids to identical molecules.
#[derive(Clone, Debug, Default)]
pub struct SpeciesRegistry {
    molecules: Vec<Molecule>,
    canon: Vec<CanonicalForm>,
}

/// A dense handle to a molecular species. What cells actually store.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub struct SpeciesId(pub u32);

impl SpeciesRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    /// Look up a molecule, registering it if new. The single entry point by which
    /// the world learns that a substance exists.
    pub fn intern(&mut self, mol: Molecule) -> SpeciesId {
        let c = mol.canonical();
        if let Some(i) = self.canon.iter().position(|x| *x == c) {
            return SpeciesId(i as u32);
        }
        let id = self.molecules.len() as u32;
        self.molecules.push(mol);
        self.canon.push(c);
        SpeciesId(id)
    }

    /// Look up without registering. `None` if this molecule has never been seen.
    pub fn lookup(&self, mol: &Molecule) -> Option<SpeciesId> {
        let c = mol.canonical();
        self.canon.iter().position(|x| *x == c).map(|i| SpeciesId(i as u32))
    }

    pub fn get(&self, id: SpeciesId) -> Option<&Molecule> {
        self.molecules.get(id.0 as usize)
    }

    pub fn canonical(&self, id: SpeciesId) -> Option<&CanonicalForm> {
        self.canon.get(id.0 as usize)
    }

    pub fn len(&self) -> usize {
        self.molecules.len()
    }
    pub fn is_empty(&self) -> bool {
        self.molecules.is_empty()
    }
    pub fn iter(&self) -> impl Iterator<Item = (SpeciesId, &Molecule)> {
        self.molecules.iter().enumerate().map(|(i, m)| (SpeciesId(i as u32), m))
    }

    /// Serialize the registry — every molecule, in id order.
    ///
    /// Under Phase 3's declared chemistry this was unnecessary: the registry was
    /// a pure function of config, rebuilt identically at every boot. **Discovery
    /// changes that.** Once a reaction can intern a molecule no config ever
    /// named, the mapping from id to graph is *history* — a fact about what this
    /// particular world has done — and history must ride in the snapshot, or a
    /// resumed world would stare at populations of species it cannot name.
    ///
    /// Wire format, all little-endian: `u32` molecule count; then per molecule
    /// `u16` atom count, one `u8` atomic number per atom, `u16` bond count, and
    /// `(u16 a, u16 b, u8 order)` per bond. Atom order within a molecule is
    /// preserved exactly — ids must survive the round trip, and re-canonicalising
    /// on load must reproduce byte-identical canonical forms, which a test pins.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(&(self.molecules.len() as u32).to_le_bytes());
        for m in &self.molecules {
            out.extend_from_slice(&(m.atom_count() as u16).to_le_bytes());
            for a in m.atoms() {
                out.push(a.z.0);
            }
            out.extend_from_slice(&(m.bonds().len() as u16).to_le_bytes());
            for b in m.bonds() {
                out.extend_from_slice(&b.a.to_le_bytes());
                out.extend_from_slice(&b.b.to_le_bytes());
                out.push(b.order.slots());
            }
        }
        out
    }

    /// Rebuild a registry from [`SpeciesRegistry::to_bytes`] output. Errors on any
    /// malformation — a truncated registry means a snapshot that cannot be
    /// trusted, and limping past it would let cells reference phantom species.
    pub fn from_bytes(bytes: &[u8]) -> Result<SpeciesRegistry, String> {
        use crate::molecule::{Atom, Bond, BondOrder};
        let mut p = 0usize;
        let take = |p: &mut usize, n: usize| -> Result<&[u8], String> {
            if *p + n > bytes.len() {
                return Err("species registry truncated".into());
            }
            let s = &bytes[*p..*p + n];
            *p += n;
            Ok(s)
        };
        let count = u32::from_le_bytes(take(&mut p, 4)?.try_into().unwrap()) as usize;
        let mut reg = SpeciesRegistry::new();
        for i in 0..count {
            let na = u16::from_le_bytes(take(&mut p, 2)?.try_into().unwrap()) as usize;
            let mut atoms = Vec::with_capacity(na);
            for _ in 0..na {
                atoms.push(Atom { z: crate::element::Z(take(&mut p, 1)?[0]) });
            }
            let nb = u16::from_le_bytes(take(&mut p, 2)?.try_into().unwrap()) as usize;
            let mut bonds = Vec::with_capacity(nb);
            for _ in 0..nb {
                let a = u16::from_le_bytes(take(&mut p, 2)?.try_into().unwrap());
                let b = u16::from_le_bytes(take(&mut p, 2)?.try_into().unwrap());
                let o = take(&mut p, 1)?[0];
                let order = BondOrder::from_u8(o)
                    .ok_or_else(|| format!("molecule {}: bond order {} invalid", i, o))?;
                bonds.push(Bond::new(a, b, order));
            }
            let id = reg.intern(Molecule::new(atoms, bonds));
            if id.0 as usize != i {
                return Err(format!(
                    "molecule {} deduplicated to id {} on load — the registry contained \
                     an isomorphic pair, which discovery can never produce",
                    i, id.0
                ));
            }
        }
        if p != bytes.len() {
            return Err("species registry has trailing bytes".into());
        }
        Ok(reg)
    }
}

/// One side of a reaction: a species and how many of it participate.
#[derive(Clone, Copy, Debug)]
pub struct Term {
    pub species: SpeciesId,
    pub count: u32,
}

/// A reaction rule.
///
/// Reactants become products. The activation energy gates the rate (how high the
/// hill is before the reaction can proceed); the enthalpy is how much energy is
/// released (negative, exothermic) or absorbed (positive, endothermic). Enthalpy
/// is *derived* from the bond energies of reactants and products at construction —
/// it is not a free parameter, because if it were, one could hand-author a reaction
/// that violated energy conservation, which is the whole thing we refuse to permit.
#[derive(Clone, Debug)]
pub struct Reaction {
    pub reactants: Vec<Term>,
    pub products: Vec<Term>,

    /// Activation energy, J·mol⁻¹. The height of the kinetic barrier. Governs *how
    /// fast*, never *whether* — a reaction with a huge barrier is slow, not
    /// forbidden. Always ≥ 0, and ≥ the endothermicity (you cannot have a barrier
    /// lower than the net uphill climb).
    pub activation_j: f64,

    /// Enthalpy of reaction, J per reaction event, as exact energy. Negative =
    /// exothermic (heats the cell). This is the quantity that couples chemistry to
    /// the Phase 2 thermal field.
    ///
    /// **Dynamic-range note.** `Energy` is quantised to the microjoule, and a
    /// single molecular reaction releases *far* less than a microjoule — the
    /// enthalpy of one H–O bond forming is ~5×10⁻¹⁹ J. Stored here alone it would
    /// round to zero, and a whole planet's combustion would release no heat. So the
    /// reactor does not use this field for the actual energy transfer; it uses
    /// [`Reaction::enthalpy_j`], the full-precision per-event value in joules, and
    /// quantises only the *aggregate* `events × enthalpy_j` once. This field is kept
    /// as a convenience mirror for code that wants an `Energy` and is working at
    /// scales where the quantum is irrelevant.
    ///
    /// **In practice, for any bond-derived reaction this field is zero.** A
    /// correct per-event enthalpy is a ten-millionth of the quantum. Read
    /// `enthalpy_j` and mean it; anything comparing this field against zero is
    /// asking a question it cannot answer.
    pub enthalpy: Energy,

    /// Enthalpy per event, in joules, at full `f64` precision.
    ///
    /// The reactor multiplies this by the (large) event count and quantises the
    /// product, so sub-microjoule per-event energies still sum to exact,
    /// non-zero heat. This is the honest way to couple molecule-scale chemistry to
    /// a microjoule-quantised thermal field: quantise the sum, never the summand.
    pub enthalpy_j: f64,

    /// Pre-exponential factor (the "A" in Arrhenius), s⁻¹ or concentration-adjusted.
    /// Sets the overall timescale; folds in collision frequency and orientation.
    /// The **structural** part of the Arrhenius prefactor. Its units depend on
    /// the reaction's molecularity, which is why it must never be read directly
    /// — call [`Reaction::rate_coefficient`], which is the one place that knows
    /// how to reconstitute `A(T)`:
    ///
    ///   * unimolecular — a dimensionless steric factor; `A = this · k_B·T/h`
    ///   * bimolecular  — m³·molecule⁻¹·s⁻¹·K^(−½); `A = this · √T`
    ///   * termolecular — the same, times an encounter volume
    ///
    /// See [`crate::kinetics`] for why a single number could not serve all
    /// three, and what was wrong before it did.
    pub pre_exponential: f64,
}

impl Reaction {
    /// The Arrhenius rate coefficient at a given temperature:
    ///
    /// ```text
    /// k(T) = A · exp(−Eₐ / RT)
    /// ```
    ///
    /// This single expression is why heating speeds reactions up — and why it does
    /// so *nonlinearly and dramatically*, roughly doubling the rate every 10 K near
    /// room temperature. It is the coupling that makes the chemistry–heat loop
    /// lively rather than languid: a little warming can wake a dormant reaction,
    /// which releases heat, which wakes more. That positive feedback is the seed of
    /// everything from combustion to metabolism.
    #[inline]
    pub fn rate_coefficient(&self, temperature_k: f64) -> f64 {
        if temperature_k <= 0.0 {
            return 0.0;
        }
        let a = crate::kinetics::prefactor_at(
            self.pre_exponential,
            self.molecularity(),
            temperature_k,
        );
        a * (-self.activation_j / (R_GAS * temperature_k)).exp()
    }

    /// How many molecules have to meet for this reaction to happen — the sum of
    /// the reactant stoichiometry. One is a molecule falling apart on its own;
    /// two is a collision; three is a collision that also has to find a third
    /// body, which is why it is rare.
    ///
    /// Derived rather than stored, so it can never disagree with the terms.
    pub fn molecularity(&self) -> u32 {
        self.reactants.iter().map(|t| t.count).sum()
    }

    /// Net atom change of this reaction, indexed by element. Must be all-zero, or
    /// the reaction creates or destroys atoms. Checked at construction.
    pub fn atom_balance(
        &self,
        table: &PeriodicTable,
        reg: &SpeciesRegistry,
    ) -> Vec<i128> {
        let mut bal = vec![0i128; table.len()];
        for t in &self.reactants {
            if let Some(m) = reg.get(t.species) {
                for (i, c) in m.formula(table).into_iter().enumerate() {
                    bal[i] -= c * t.count as i128;
                }
            }
        }
        for t in &self.products {
            if let Some(m) = reg.get(t.species) {
                for (i, c) in m.formula(table).into_iter().enumerate() {
                    bal[i] += c * t.count as i128;
                }
            }
        }
        bal
    }

    /// Is this reaction atom-balanced? A reaction that fails this must never enter
    /// the simulation — it is a hole in conservation waiting to be exploited.
    pub fn is_atom_balanced(&self, table: &PeriodicTable, reg: &SpeciesRegistry) -> bool {
        self.atom_balance(table, reg).iter().all(|&b| b == 0)
    }
}

/// A model that assigns an energy to a bond, and thereby an enthalpy to a reaction.
///
/// Bond energies are *derived*, not tabulated pair-by-pair, from the elements'
/// electronegativities and the bond order. This is a deliberate, heavy
/// simplification of real quantum chemistry — a real bond energy depends on the
/// whole molecular environment — but it captures the essential trend (bonds between
/// atoms of very different electronegativity, like O–H, are strong; the energy rises
/// with bond order) with a handful of parameters rather than a giant table. It is
/// the same bet the whole crate makes: preserve the causal trend, pay almost
/// nothing, and let Phase 5 tell us if a refinement is needed.
#[derive(Clone, Debug)]
pub struct BondEnergyModel {
    /// Base bond energy per shared electron pair, J·mol⁻¹.
    pub base_per_order: f64,
    /// How much an electronegativity difference strengthens a bond, J·mol⁻¹ per
    /// (Pauling unit)². (Pauling's own relation is quadratic in the difference.)
    pub electroneg_coefficient: f64,
}

impl Default for BondEnergyModel {
    fn default() -> Self {
        // Values in the right ballpark for real single-bond energies (a C–C single
        // bond is ~350 kJ/mol, an O–H ~460). Tuned, in the validation fixture, so
        // the engine reproduces a known equilibrium; not physically fundamental.
        BondEnergyModel { base_per_order: 200_000.0, electroneg_coefficient: 100_000.0 }
    }
}

impl BondEnergyModel {
    /// The **potential energy** of a molecule, relative to its free atoms.
    ///
    /// # The sign convention, and why the first version was backwards
    ///
    /// A bond *releases* energy when it forms — that is what "stable" means. So a
    /// molecule with many strong bonds sits *lower* on the energy scale than the
    /// same atoms floating free, not higher. My first attempt tracked "energy
    /// stored in a bond" as a positive quantity and made strongly-bonded molecules
    /// *higher* in energy, which inverted every enthalpy: forming a strong O–H bond
    /// came out endothermic, the exact opposite of reality (hydrogen and oxygen do
    /// not gently absorb heat to make water).
    ///
    /// The fix is to measure potential energy from the free-atom baseline of zero
    /// and count each bond as a *negative* contribution — the deeper the well, the
    /// more stable the molecule. A reaction is then exothermic exactly when its
    /// products sit lower than its reactants, which is the physical definition.
    ///
    /// This is dimensional bookkeeping done right, not a new model: the same bond
    /// strengths, entered with the sign that makes "stable = low energy" true.
    pub fn molecule_energy(&self, mol: &Molecule, table: &PeriodicTable) -> f64 {
        let mut binding = 0.0;
        for b in mol.bonds() {
            let za = mol.atoms()[b.a as usize].z;
            let zb = mol.atoms()[b.b as usize].z;
            let (Some(ea), Some(eb)) = (table.get(za), table.get(zb)) else { continue };
            let d = (ea.electronegativity() - eb.electronegativity()).abs();
            // Strength of this bond (always positive): base, plus a polarity bonus.
            let strength = self.base_per_order + self.electroneg_coefficient * d * d;
            binding += strength * b.order.slots() as f64;
        }
        // Potential energy is the *negative* of total binding: strongly-bonded
        // molecules are deep in the well.
        -binding
    }

    /// Enthalpy of a reaction, **per mole**: energy(products) − energy(reactants).
    ///
    /// Negative ⇒ products sit lower than reactants ⇒ energy released ⇒ exothermic.
    /// Deriving it as a difference of potential energies makes conservation
    /// automatic: the heat delivered to the thermal field is, by construction,
    /// exactly the depth the system dropped in its own potential well.
    ///
    /// # Per mole, and it says so now
    ///
    /// `base_per_order` is quoted in J·mol⁻¹ because that is how bond energies
    /// are tabulated, so everything built from it is molar. This function used to
    /// promise "per event" in its documentation while returning a molar `Energy`,
    /// and two of its three callers believed the prose: `with_derived_enthalpy`
    /// and the declared-reaction parser both assigned it straight to
    /// `enthalpy_j`, spending a mole of bond energy on a single molecule — every
    /// declared reaction without an explicit `enthalpy_j` delivered 6.022e23
    /// times too much heat. `network.rs` was right the whole time; it named its
    /// variable `dh_molar` and divided.
    ///
    /// The unit now lives in the name and in the return type — a bare `f64` of
    /// joules per mole rather than an `Energy`, because `Energy` is quantised
    /// attojoules and a correct per-event enthalpy (~1e-19 J) is a *fraction* of
    /// one attojoule. Converting first and quantising second is the only order
    /// that survives; `Reaction::enthalpy_j` exists for exactly this reason.
    ///
    /// Callers must divide by [`crate::network::AVOGADRO`].
    pub fn reaction_enthalpy_molar(
        &self,
        reactants: &[Term],
        products: &[Term],
        reg: &SpeciesRegistry,
        table: &PeriodicTable,
    ) -> f64 {
        let side = |terms: &[Term]| -> f64 {
            terms
                .iter()
                .filter_map(|t| reg.get(t.species).map(|m| self.molecule_energy(m, table) * t.count as f64))
                .sum()
        };
        let reactant_e = side(reactants);
        let product_e = side(products);
        // ΔH = E(products) − E(reactants), joules per mole. Lower products ⇒
        // negative ⇒ exothermic.
        //
        // No division here on purpose. Each caller knows whether it wants molar
        // or per-event, and making them ask for it is what stops the next reader
        // from guessing — the previous version guessed, in prose, and was wrong
        // for six phases while energy conservation held perfectly throughout.
        // A conservation check cannot see a dimensional error: the heat was
        // exactly accounted for on both sides and simply the wrong size.
        // Invariants catch bookkeeping mistakes; only an external anchor catches
        // a units mistake.
        product_e - reactant_e
    }
}
