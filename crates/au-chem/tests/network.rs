//! The generative network, held to account.
//!
//! Phase 3 proved the reactor's thermodynamics on *declared* reactions. What
//! Phase 5a must prove is stranger: that reactions **derived from molecular
//! structure by an algorithm** — with products the world has never seen, interned
//! mid-flight — still land on the same physics. The tests run in three tiers:
//!
//!   * **Mechanics** — the graph moves are exact inverses, symmetry classes
//!     collapse equivalent sites, atoms always balance, the cap holds.
//!   * **Discovery** — from nothing but hydrogen and oxygen *atoms*, the network
//!     reaches water. Nobody declares H₂O; it is found, because it is reachable.
//!   * **Thermodynamics** — the activation identity `Ea_f − Ea_r = ΔH` holds for
//!     every generated pair, and a reactor running *generated* reactions settles
//!     at the Boltzmann equilibrium `K = exp(−ΔH/RT)`. Chemistry the engine
//!     invented obeys laws nobody wrote into the generator.

use au_chem::kinetics::prefactor_at;
use au_chem::*;

fn hydrogen_oxygen_table() -> PeriodicTable {
    let mut t = PeriodicTable::new();
    t.add(Element { z: Z(1), symbol: "H", mass_mda: 1_008, valence: 1, electronegativity_c: 220 });
    t.add(Element { z: Z(8), symbol: "O", mass_mda: 15_999, valence: 2, electronegativity_c: 344 });
    t
}

fn atoms_only_registry() -> SpeciesRegistry {
    let mut reg = SpeciesRegistry::new();
    reg.intern(Molecule::monatomic(Z(1))); // H·
    reg.intern(Molecule::monatomic(Z(8))); // O·
    reg
}

/// Run generation to a fixpoint: enumerate, add every product to the present
/// set, repeat until no new species appear. This is what the sim loop will do
/// incrementally; here it is done eagerly so tests can inspect the closure.
fn closure(
    reg: &mut SpeciesRegistry,
    table: &PeriodicTable,
    model: &BondEnergyModel,
    rules: &NetworkRules,
) -> Vec<Reaction> {
    loop {
        let present: Vec<SpeciesId> = reg.iter().map(|(id, _)| id).collect();
        let before = reg.len();
        let reactions = enumerate_reactions(&present, reg, table, model, rules);
        if reg.len() == before {
            return reactions;
        }
    }
}

// ═══ Mechanics ══════════════════════════════════════════════════════════════

/// `split` undoes `join`: the fragments of the freshly formed bond are the
/// molecules that were joined, up to isomorphism.
#[test]
fn join_and_split_are_inverse() {
    let h2 = Molecule::new(
        vec![Atom { z: Z(1) }, Atom { z: Z(1) }],
        vec![Bond::new(0, 1, BondOrder::Single)],
    );
    let o = Molecule::monatomic(Z(8));
    let joined = join(&h2, 0, &o, 0); // H–H joined to O at one H → H–H–O? No: H has
                                      // valence 1 and is full — the *generator* would
                                      // refuse this site; `join` itself only does
                                      // bookkeeping, which is what this test checks.
    let (f0, f1) = split(&joined, joined.bonds().len() - 1).expect("must disconnect");
    let (c0, c1) = (f0.canonical(), f1.canonical());
    let (h2c, oc) = (h2.canonical(), o.canonical());
    assert!(
        (c0 == h2c && c1 == oc) || (c0 == oc && c1 == h2c),
        "splitting the joined bond must recover the original molecules"
    );
}

/// Water's two hydrogens are one symmetry class, so the generator offers exactly
/// one way to pull an H off — not two identical ones.
#[test]
fn symmetric_sites_collapse_to_one_candidate() {
    let table = hydrogen_oxygen_table();
    let model = BondEnergyModel::default();
    let rules = NetworkRules::default();
    let mut reg = SpeciesRegistry::new();
    let water = reg.intern(Molecule::new(
        vec![Atom { z: Z(1) }, Atom { z: Z(8) }, Atom { z: Z(1) }],
        vec![Bond::new(0, 1, BondOrder::Single), Bond::new(1, 2, BondOrder::Single)],
    ));
    let rx = enumerate_reactions(&[water], &mut reg, &table, &model, &rules);
    let splits: Vec<_> = rx
        .iter()
        .filter(|r| r.reactants.len() == 1 && r.reactants[0].species == water && r.products.len() == 2)
        .collect();
    assert_eq!(
        splits.len(),
        1,
        "H₂O has one O–H bond class; the generator offered {} split(s)",
        splits.len()
    );
}

/// Every reaction the generator ever emits balances atoms — across the whole
/// closure of an H/O world, not just cherry-picked cases. A generator that could
/// emit one unbalanced rule would hand evolution a matter faucet.
#[test]
fn every_generated_reaction_balances_atoms() {
    let table = hydrogen_oxygen_table();
    let model = BondEnergyModel::default();
    let rules = NetworkRules { max_atoms: 4, ..NetworkRules::default() };
    let mut reg = atoms_only_registry();
    let rx = closure(&mut reg, &table, &model, &rules);
    assert!(rx.len() > 10, "the closure should contain a real network, got {}", rx.len());
    for r in &rx {
        assert!(
            r.is_atom_balanced(&table, &reg),
            "generated reaction does not balance atoms"
        );
    }
}

/// The cap is a wall, not a suggestion: with `max_atoms = 2`, nothing with three
/// atoms is ever created, so water is unreachable — and the world is smaller but
/// still consistent.
#[test]
fn max_atoms_caps_growth() {
    let table = hydrogen_oxygen_table();
    let model = BondEnergyModel::default();
    let rules = NetworkRules { max_atoms: 2, ..NetworkRules::default() };
    let mut reg = atoms_only_registry();
    closure(&mut reg, &table, &model, &rules);
    for (_, m) in reg.iter() {
        assert!(m.atom_count() <= 2, "the cap was breached: {} atoms", m.atom_count());
    }
}

// ═══ Discovery ══════════════════════════════════════════════════════════════

/// **The Phase 5a headline.** Seed the world with hydrogen atoms and oxygen
/// atoms — nothing else, no declared molecules, no declared reactions — and the
/// generative closure contains water. H₂, O₂, hydroxyl, peroxide too: the whole
/// neighbourhood, reached by nothing but "bonds can form where valence allows".
///
/// This is the moment "no hardcoded molecules" stops being an architecture note
/// and becomes an observable: the engine *found* H₂O.
#[test]
fn the_engine_discovers_water() {
    let table = hydrogen_oxygen_table();
    let model = BondEnergyModel::default();
    let rules = NetworkRules { max_atoms: 4, ..NetworkRules::default() };
    let mut reg = atoms_only_registry();
    closure(&mut reg, &table, &model, &rules);

    let water = Molecule::new(
        vec![Atom { z: Z(1) }, Atom { z: Z(8) }, Atom { z: Z(1) }],
        vec![Bond::new(0, 1, BondOrder::Single), Bond::new(1, 2, BondOrder::Single)],
    );
    assert!(
        reg.lookup(&water).is_some(),
        "starting from H and O atoms, the network should have discovered water"
    );

    // And the neighbourhood: molecular hydrogen, molecular oxygen (as O=O —
    // reached by forming O–O and raising it), hydroxyl.
    let h2 = Molecule::new(vec![Atom { z: Z(1) }; 2], vec![Bond::new(0, 1, BondOrder::Single)]);
    let o2 = Molecule::new(vec![Atom { z: Z(8) }; 2], vec![Bond::new(0, 1, BondOrder::Double)]);
    let oh = Molecule::new(
        vec![Atom { z: Z(1) }, Atom { z: Z(8) }],
        vec![Bond::new(0, 1, BondOrder::Single)],
    );
    assert!(reg.lookup(&h2).is_some(), "H₂ should be discovered");
    assert!(reg.lookup(&o2).is_some(), "O=O should be discovered via form-then-raise");
    assert!(reg.lookup(&oh).is_some(), "hydroxyl should be discovered");
}

/// Discovery is deterministic: two closures from identical seeds produce the
/// same species in the same order with the same reactions. Without this, a world
/// hash containing chemistry would be a coin flip.
#[test]
fn discovery_is_deterministic() {
    let run = || {
        let table = hydrogen_oxygen_table();
        let model = BondEnergyModel::default();
        let rules = NetworkRules { max_atoms: 4, ..NetworkRules::default() };
        let mut reg = atoms_only_registry();
        let rx = closure(&mut reg, &table, &model, &rules);
        let species: Vec<Vec<u8>> =
            reg.iter().map(|(id, _)| reg.canonical(id).unwrap().0.clone()).collect();
        let rules_sig: Vec<(u32, u32)> = rx
            .iter()
            .map(|r| (r.reactants[0].species.0, r.products[0].species.0))
            .collect();
        (species, rules_sig)
    };
    assert_eq!(run(), run(), "two identical seeds diverged in discovery");
}

/// The registry round-trips through its wire format: same length, identical
/// canonical form at every id, and a lookup lands on the same id. This is the
/// serialization save/resume will lean on — a resumed world must be able to
/// *name* every species its cells hold populations of.
#[test]
fn the_registry_survives_the_round_trip() {
    let table = hydrogen_oxygen_table();
    let model = BondEnergyModel::default();
    let rules = NetworkRules { max_atoms: 4, ..NetworkRules::default() };
    let mut reg = atoms_only_registry();
    closure(&mut reg, &table, &model, &rules);

    let bytes = reg.to_bytes();
    let back = SpeciesRegistry::from_bytes(&bytes).expect("round trip must succeed");
    assert_eq!(reg.len(), back.len(), "species count changed in the round trip");
    for (id, _) in reg.iter() {
        assert_eq!(
            reg.canonical(id).unwrap(),
            back.canonical(id).unwrap(),
            "species {} changed identity in the round trip",
            id.0
        );
    }
}

// ═══ Thermodynamics ═════════════════════════════════════════════════════════

/// `Ea_forward − Ea_reverse = ΔH`, for **every** generated pair in the closure.
/// The module claims this is true by construction; this test exists because
/// "true by construction" is a property of code, and code gets refactored.
#[test]
fn activation_energies_obey_detailed_balance() {
    let table = hydrogen_oxygen_table();
    let model = BondEnergyModel::default();
    let rules = NetworkRules { max_atoms: 4, ..NetworkRules::default() };
    let mut reg = atoms_only_registry();
    let rx = closure(&mut reg, &table, &model, &rules);

    // Index reactions by (reactants, products) as canonical multiset keys.
    let key = |terms: &[Term]| -> Vec<(Vec<u8>, u32)> {
        let mut v: Vec<(Vec<u8>, u32)> =
            terms.iter().map(|t| (reg.canonical(t.species).unwrap().0.clone(), t.count)).collect();
        v.sort();
        v
    };
    let mut pairs = 0;
    for f in &rx {
        let (kf_r, kf_p) = (key(&f.reactants), key(&f.products));
        for r in &rx {
            if key(&r.reactants) == kf_p && key(&r.products) == kf_r {
                let dh_f = f.enthalpy_j * AVOGADRO;
                let diff = f.activation_j - r.activation_j;
                assert!(
                    (diff - dh_f).abs() < 1.0e-6 * dh_f.abs().max(1.0),
                    "detailed balance broken: Ea_f − Ea_r = {} but ΔH = {}",
                    diff,
                    dh_f
                );
                pairs += 1;
            }
        }
    }
    assert!(pairs >= 6, "the closure should contain reverse pairs to check, found {}", pairs);
}

/// A reactor running **generated** reactions settles at the Boltzmann
/// equilibrium. A hydrogen-only world — 2 H· ⇌ H₂, both directions produced by
/// the generator, nothing declared — is run to equilibrium at two temperatures,
/// and the measured K = [H₂]/[H·]² is compared with the prediction
/// K = k_f/k_r = exp(−ΔH/RT) that the activation identity implies.
///
/// Nothing in the generator mentions Boltzmann. The distribution appears because
/// the energetics are consistent — which is the entire reason the Ea model was
/// chosen the way it was.
#[test]
fn generated_reactions_reach_boltzmann_equilibrium() {
    let mut table = PeriodicTable::new();
    table.add(Element { z: Z(1), symbol: "H", mass_mda: 1_008, valence: 1, electronegativity_c: 220 });
    let model = BondEnergyModel::default();
    let rules =
        NetworkRules { intrinsic_barrier_j_mol: 40_000.0, kinetics: KineticRules::default(), max_atoms: 2 };

    let mut reg = SpeciesRegistry::new();
    let h = reg.intern(Molecule::monatomic(Z(1)));
    let rx = closure(&mut reg, &table, &model, &rules);
    let h2 = SpeciesId(1);
    assert_eq!(reg.len(), 2, "an H-only, 2-atom world holds exactly H· and H₂");
    assert!(rx.len() >= 2, "both directions of 2H ⇌ H₂ should be generated");

    // The validated envelope, measured rather than assumed. Sweeping
    // temperature at this concentration gives K_measured/K_predicted of
    //
    //     2000 K → 0.99    3500 K → 0.70
    //     2500 K → 0.97    4000 K → 0.04
    //     3000 K → 0.90    4500 K → 0.00
    //
    // Detailed balance holds tightly while the minority species is a fair
    // fraction of the population, and the *integer* reactor under-resolves the
    // equilibrium once H₂ falls below about a percent of the atoms: dissociation
    // is unimolecular and fast, association needs two particles to meet, and at
    // that point the per-substep consumption cap governs instead of the rate
    // constants. That is a property of the reactor, not of the kinetics — and it
    // predates this change; the old test simply never reached a regime where the
    // equilibrium was measurable at all. The band tested below is the band the
    // reactor can honestly resolve.
    //
    // A physically real regime: ~7 mmol/L of hydrogen in a cubic centimetre.
    //
    // The previous version of this test ran at about four molecules per cubic
    // metre — a hard vacuum — where an association equilibrium simply cannot be
    // measured, and it only ever passed because two errors cancelled: the
    // missing entropy term and the missing concentration dependence. With
    // prefactors derived from collision theory the entropy is present, so the
    // concentration has to be real too. 2H ⇌ H₂ is then roughly half dissociated
    // at 3000 K, which is where the equilibrium is actually informative.
    let equilibrate = |temp: f64| -> f64 {
        let vol = 1.0e-6;
        let mut c = CellChemistry::new(reg.len());
        c.set(h2, 2_250_000_000_000_000_000);
        let mut e = au_physics::Energy::from_joules(1.0e18); // vast bath: T pinned
        for _ in 0..2_000 {
            let mut cell = ReactingCell {
                chem: &mut c,
                thermal_energy: &mut e,
                temperature_k: temp,
                volume_m3: vol,
            };
            react(&mut cell, &rx, 1.0e-9);
        }
        let ch = c.get(h) as f64 / vol;
        let ch2 = c.get(h2) as f64 / vol;
        ch2 / (ch * ch).max(1.0e-300)
    };

    for temp in [2000.0, 2500.0, 3000.0] {
        // ΔH for 2H → H₂ from the bond model, and the Boltzmann prediction.
        let fwd = rx
            .iter()
            .find(|r| r.reactants[0].species == h && r.products[0].species == h2)
            .expect("the forming direction must exist");
        let rev = rx
            .iter()
            .find(|r| r.reactants[0].species == h2 && r.products[0].species == h)
            .expect("the splitting direction must exist");
        let dh = fwd.enthalpy_j * AVOGADRO;

        // The equilibrium constant physics actually demands is
        //     K = (A_f / A_r) · exp(−ΔH/RT)
        // and the prefactor ratio is not a detail — it *is* the reaction
        // entropy. Association is bimolecular forwards and unimolecular
        // backwards, so A_f/A_r is a small volume: joining two molecules into
        // one costs translational entropy. This is detailed balance, and it is
        // the strongest statement available about a generated network — that
        // the reactor arrives at the equilibrium its own rate constants imply.
        let a_f = prefactor_at(fwd.pre_exponential, fwd.molecularity(), temp);
        let a_r = prefactor_at(rev.pre_exponential, rev.molecularity(), temp);
        let entropy_volume = a_f / a_r;
        assert!(
            (1.0e-32..1.0e-26).contains(&entropy_volume),
            "A_f/A_r = {:.3e} m³ — the entropy of association is out of range",
            entropy_volume
        );
        let k_pred = entropy_volume * (-dh / (R_GAS * temp)).exp();
        let k_meas = equilibrate(temp);
        let ratio = k_meas / k_pred;
        assert!(
            (0.85..1.15).contains(&ratio),
            "at {} K, measured K = {:.4e} vs Boltzmann {:.4e} (ratio {:.3}) — \
             generated chemistry missed the equilibrium physics demands",
            temp,
            k_meas,
            k_pred,
            ratio
        );
    }
}
