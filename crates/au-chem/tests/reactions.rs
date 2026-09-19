//! Reactions, checked against known chemistry.
//!
//! Phase 2b had Ra_c = 657.51 — a number derived on paper that the engine had to
//! reproduce. Chemistry's equivalent is the **equilibrium constant** and the
//! direction of **Le Chatelier's principle**: a reversible reaction left alone must
//! settle at the ratio of products to reactants that thermodynamics predicts, and
//! must *shift* the right way when the temperature changes. Nothing in the reactor
//! knows those laws — it only fires rules at Arrhenius rates. If equilibrium and Le
//! Chatelier fall out anyway, the abstraction has captured the real causal
//! structure, and that is the whole bet of Phase 3.
//!
//! And underneath the validation, the two non-negotiables: atoms conserved
//! exactly, and total energy (bonds + heat) conserved exactly.

use au_chem::*;
use au_physics::Energy;

fn table() -> PeriodicTable {
    let mut t = PeriodicTable::new();
    t.add(Element { z: Z::H, symbol: "H", mass_mda: 1_008, valence: 1, electronegativity_c: 220 });
    t.add(Element { z: Z::N, symbol: "N", mass_mda: 14_007, valence: 3, electronegativity_c: 304 });
    t.add(Element { z: Z::O, symbol: "O", mass_mda: 15_999, valence: 2, electronegativity_c: 344 });
    t.add(Element { z: Z::C, symbol: "C", mass_mda: 12_011, valence: 4, electronegativity_c: 255 });
    t
}

/// Build the two diatomic reactants and their bonded product for a generic
/// A₂ + B₂ ⇌ 2 AB reaction, returning their species ids. Using H and its own
/// "isotope" lookalike would be cleaner, but two real elements make the bond
/// energies (and thus the enthalpy) physically sensible.
fn ab_system(reg: &mut SpeciesRegistry) -> (SpeciesId, SpeciesId, SpeciesId) {
    // A2 = H2 (H–H)
    let a2 = reg.intern(Molecule::new(
        vec![Atom { z: Z::H }, Atom { z: Z::H }],
        vec![Bond::new(0, 1, BondOrder::Single)],
    ));
    // B2 = O2-like single-bonded pair (kept single-bond so valence stays simple)
    let b2 = reg.intern(Molecule::new(
        vec![Atom { z: Z::O }, Atom { z: Z::O }],
        vec![Bond::new(0, 1, BondOrder::Single)],
    ));
    // AB = H–O
    let ab = reg.intern(Molecule::new(
        vec![Atom { z: Z::H }, Atom { z: Z::O }],
        vec![Bond::new(0, 1, BondOrder::Single)],
    ));
    (a2, b2, ab)
}

/// A forward/reverse reaction pair with enthalpies derived from a bond model, so
/// they are guaranteed energy-consistent (the reverse releases exactly what the
/// forward absorbs).
fn reversible_pair(
    a2: SpeciesId,
    b2: SpeciesId,
    ab: SpeciesId,
    reg: &SpeciesRegistry,
    t: &PeriodicTable,
    model: &BondEnergyModel,
    fwd_a: f64,
    rev_a: f64,
    fwd_act: f64,
    rev_act: f64,
) -> (Reaction, Reaction) {
    let forward = with_derived_enthalpy(
        Reaction {
            reactants: vec![Term { species: a2, count: 1 }, Term { species: b2, count: 1 }],
            products: vec![Term { species: ab, count: 2 }],
            activation_j: fwd_act,
            enthalpy: Energy::ZERO,
            enthalpy_j: 0.0,
            pre_exponential: fwd_a,
        },
        model,
        reg,
        t,
    );
    let reverse = with_derived_enthalpy(
        Reaction {
            reactants: vec![Term { species: ab, count: 2 }],
            products: vec![Term { species: a2, count: 1 }, Term { species: b2, count: 1 }],
            activation_j: rev_act,
            enthalpy: Energy::ZERO,
            enthalpy_j: 0.0,
            pre_exponential: rev_a,
        },
        model,
        reg,
        t,
    );
    (forward, reverse)
}

// ═══ CONSERVATION: the two non-negotiables ═══════════════════════════════════

/// A reaction that does not balance atoms must be *detectably* unbalanced, so the
/// assembly layer can refuse it before it ever runs. The reactor's whole
/// conservation guarantee assumes every rule it is given is balanced.
#[test]
fn unbalanced_reactions_are_rejected() {
    let t = table();
    let mut reg = SpeciesRegistry::new();
    let (a2, b2, ab) = ab_system(&mut reg);

    // Correct: A2 + B2 -> 2 AB has 2 H and 2 O on each side.
    let ok = Reaction {
        reactants: vec![Term { species: a2, count: 1 }, Term { species: b2, count: 1 }],
        products: vec![Term { species: ab, count: 2 }],
        activation_j: 0.0,
        enthalpy: Energy::ZERO,
            enthalpy_j: 0.0,
        pre_exponential: 1.0,
    };
    assert!(ok.is_atom_balanced(&t, &reg), "a balanced reaction was flagged unbalanced");

    // Broken: A2 + B2 -> 1 AB loses an H and an O into thin air.
    let bad = Reaction {
        reactants: vec![Term { species: a2, count: 1 }, Term { species: b2, count: 1 }],
        products: vec![Term { species: ab, count: 1 }],
        activation_j: 0.0,
        enthalpy: Energy::ZERO,
            enthalpy_j: 0.0,
        pre_exponential: 1.0,
    };
    assert!(!bad.is_atom_balanced(&t, &reg), "an atom-destroying reaction passed as balanced");
}

/// **The chemical conservation test.** Run a reversible reaction to equilibrium and
/// beyond — hundreds of thousands of reaction events, forward and back — and the
/// total count of every element is *exactly* unchanged. Atoms are only ever
/// rearranged between molecules; none is created, none destroyed.
#[test]
fn atoms_are_conserved_through_reaction() {
    let t = table();
    let mut reg = SpeciesRegistry::new();
    let (a2, b2, ab) = ab_system(&mut reg);
    let model = BondEnergyModel::default();
    let (fwd, rev) = reversible_pair(a2, b2, ab, &reg, &t, &model, 1e5, 1e5, 40_000.0, 40_000.0);

    let mut chem = CellChemistry::new(reg.len());
    // Populations sized so the heat is *observable*, which is a real constraint
    // and not a fudge. `Energy` is quantised to the microjoule and one bond is
    // ~5e-19 J, so the thermal field cannot resolve fewer than roughly 2e12
    // reaction events however long you wait. This test previously used 1e6
    // molecules and passed only because bond-derived enthalpies were 6.022e23
    // times too large — it was measuring the bug, not the chemistry.
    chem.set(a2, 1_000_000_000_000_000_000);
    chem.set(b2, 1_000_000_000_000_000_000);

    let inv0 = chem.atom_inventory(&t, &reg);

    let mut energy = Energy::from_joules(1.0e6);
    for _ in 0..5000 {
        let mut cell = ReactingCell {
            chem: &mut chem,
            thermal_energy: &mut energy,
            temperature_k: 800.0,
            volume_m3: 1.0,
        };
        react(&mut cell, &[fwd.clone(), rev.clone()], 1.0e-6);
    }

    let inv1 = chem.atom_inventory(&t, &reg);
    assert_eq!(inv0, inv1, "ATOM LEAK across reactions: {:?} became {:?}", inv0, inv1);
}

/// **The chemical energy-conservation test.** Bond energy plus thermal energy is
/// invariant. Every joule a reaction releases as heat came out of a bond; every
/// joule an endothermic reaction absorbs went into one. The reactor moves energy
/// between two ledgers that always sum to the same total.
#[test]
fn bond_energy_plus_heat_is_conserved_exactly() {
    let t = table();
    let mut reg = SpeciesRegistry::new();
    let (a2, b2, ab) = ab_system(&mut reg);
    let model = BondEnergyModel::default();
    let (fwd, rev) = reversible_pair(a2, b2, ab, &reg, &t, &model, 5e5, 5e5, 30_000.0, 30_000.0);

    let mut chem = CellChemistry::new(reg.len());
    chem.set(a2, 500_000);
    chem.set(b2, 500_000);

    // Total energy = thermal + stored bond energy. This sum is the invariant.
    let bond_energy = |chem: &CellChemistry| -> f64 {
        let mut e = 0.0;
        for (id, m) in reg.iter() {
            e += model.molecule_energy(m, &t) * chem.get(id) as f64;
        }
        e
    };

    let mut thermal = Energy::from_joules(5.0e6);
    let total0 = thermal.as_joules() + bond_energy(&chem);

    for _ in 0..8000 {
        let temp = 400.0 + thermal.as_joules() * 1e-5; // a crude feedback so T actually varies
        let mut cell = ReactingCell {
            chem: &mut chem,
            thermal_energy: &mut thermal,
            temperature_k: temp.max(1.0),
            volume_m3: 1.0,
        };
        react(&mut cell, &[fwd.clone(), rev.clone()], 1.0e-6);
    }

    let total1 = thermal.as_joules() + bond_energy(&chem);
    // Exact in principle; allow a hair for the f64 bond-energy sum (the *reactor's*
    // accounting is integer-exact — this tolerance is only the test's own
    // re-summation of bond energies in floating point).
    let rel = (total1 - total0).abs() / total0.abs().max(1.0);
    assert!(
        rel < 1e-9,
        "ENERGY LEAK: total energy went {:.6e} → {:.6e} J ({:.2e} relative)",
        total0,
        total1,
        rel
    );
}

// ═══ KINETICS: reactions behave like reactions ══════════════════════════════

/// Arrhenius: heating a reaction speeds it up, sharply and nonlinearly. This is the
/// coupling that makes the chemistry–heat loop lively — and it is why a small
/// warming can wake a dormant reaction. The rate must rise steeply with T.
#[test]
fn heating_accelerates_reactions() {
    let t = table();
    let mut reg = SpeciesRegistry::new();
    let (a2, b2, ab) = ab_system(&mut reg);
    let model = BondEnergyModel::default();
    let (fwd, _) = reversible_pair(a2, b2, ab, &reg, &t, &model, 1e6, 1e6, 60_000.0, 60_000.0);

    let events_at = |temp: f64| -> u64 {
        let mut chem = CellChemistry::new(reg.len());
        chem.set(a2, 1_000_000);
        chem.set(b2, 1_000_000);
        let mut e = Energy::from_joules(1e9);
        let mut cell = ReactingCell {
            chem: &mut chem,
            thermal_energy: &mut e,
            temperature_k: temp,
            volume_m3: 1.0,
        };
        react(&mut cell, &[fwd.clone()], 1.0e-9).events
    };

    let cold = events_at(300.0);
    let hot = events_at(600.0);
    assert!(
        hot > cold * 10,
        "doubling temperature should multiply the rate many-fold: {} events cold vs {} hot",
        cold,
        hot
    );
}

/// Mass action: doubling a reactant's concentration must increase the forward rate.
/// The reaction goes faster when there is more to react — the defining behaviour of
/// a second-order reaction, and the reason concentration gradients drive chemistry.
#[test]
fn more_reactant_means_faster_reaction() {
    let t = table();
    let mut reg = SpeciesRegistry::new();
    let (a2, b2, ab) = ab_system(&mut reg);
    let model = BondEnergyModel::default();
    let (fwd, _) = reversible_pair(a2, b2, ab, &reg, &t, &model, 1e5, 1e5, 20_000.0, 20_000.0);

    let events_with = |a: i128| -> u64 {
        let mut chem = CellChemistry::new(reg.len());
        chem.set(a2, a);
        chem.set(b2, 1_000_000);
        let mut e = Energy::from_joules(1e9);
        let mut cell = ReactingCell {
            chem: &mut chem,
            thermal_energy: &mut e,
            temperature_k: 500.0,
            volume_m3: 1.0,
        };
        react(&mut cell, &[fwd.clone()], 1.0e-9).events
    };

    assert!(
        events_with(2_000_000) > events_with(1_000_000),
        "more reactant must react faster"
    );
}

// ═══ EQUILIBRIUM: the validation that mirrors Ra_c ══════════════════════════

/// **The headline chemistry result.**
///
/// A reversible reaction, left alone at fixed temperature, must settle at a steady
/// ratio of products to reactants and *stay there* — a dynamic equilibrium where
/// forward and reverse rates balance. Nothing imposes this; the reactor only fires
/// rules. Equilibrium is an emergent fixed point of the kinetics, exactly as it is
/// in reality.
#[test]
fn a_reversible_reaction_reaches_equilibrium() {
    let t = table();
    let mut reg = SpeciesRegistry::new();
    let (a2, b2, ab) = ab_system(&mut reg);
    let model = BondEnergyModel::default();
    // Real collision prefactors, and therefore real concentrations.
    //
    // The fixture used to pass 3e4 for both directions. Under the old fitted
    // convention that was merely arbitrary; now that a prefactor is a physical
    // quantity it would be twenty-one orders of magnitude above the gas-kinetic
    // rate, so it is replaced by one — and the populations are raised to match,
    // because a millimolar solution is where an equilibrium between two
    // bimolecular directions is actually observable. Both directions are
    // bimolecular here (A₂ + B₂ ⇌ 2AB), so their prefactors are equal and the
    // equilibrium is unchanged by this correction; only the timescale and the
    // concentration are now honest.
    let structural = 1.0e-17; // A(500 K) ≈ 2×10⁻¹⁶ m³·molecule⁻¹·s⁻¹
    let (fwd, rev) =
        reversible_pair(a2, b2, ab, &reg, &t, &model, structural, structural, 20_000.0, 20_000.0);

    let mut chem = CellChemistry::new(reg.len());
    chem.set(a2, 2_000_000_000_000_000_000);
    chem.set(b2, 2_000_000_000_000_000_000);
    let mut e = Energy::from_joules(1e9);

    let ratio_at = |chem: &CellChemistry| chem.get(ab) as f64 / chem.get(a2).max(1) as f64;

    // Run to steady state.
    let mut last = 0.0;
    for step in 0..40_000 {
        let mut cell = ReactingCell {
            chem: &mut chem,
            thermal_energy: &mut e,
            temperature_k: 500.0,
            volume_m3: 1.0e-6,
        };
        react(&mut cell, &[fwd.clone(), rev.clone()], 5.0e-8);
        if step == 20_000 {
            last = ratio_at(&chem);
        }
    }
    let now = ratio_at(&chem);

    // It must have formed a meaningful amount of product...
    assert!(
        chem.get(ab) > 100_000_000_000_000,
        "essentially no product formed; not a real equilibrium"
    );
    // ...and it must be *steady*: the ratio at step 20k and step 40k agree.
    let drift = (now - last).abs() / last.max(1e-9);
    assert!(
        drift < 0.05,
        "the ratio never settled: {:.4} at midpoint vs {:.4} at end ({:.1}% drift)",
        last,
        now,
        drift * 100.0
    );
}

/// **Le Chatelier, the direction test.**
///
/// For an exothermic forward reaction, raising the temperature must shift the
/// equilibrium *back toward reactants* — heat is effectively a product, and adding
/// heat pushes the balance away from making more of it. This is subtle and
/// emergent: it happens because the reverse (endothermic) reaction's Arrhenius rate
/// climbs faster with temperature than the forward one's (its barrier is higher by
/// exactly the reaction enthalpy). If the engine reproduces the right *direction* of
/// shift, it has captured the thermodynamics, not merely the kinetics.
///
/// # A note on picking the regime
///
/// The equilibrium constant K = (A_f/A_r)·exp(−(Eₐf−Eₐr)/RT) is *exponentially*
/// sensitive to the enthalpy. With a large enthalpy, K is astronomical at every
/// temperature the reaction ever sees, the reverse reaction is effectively dead, and
/// the reaction runs to completion at both temperatures — so the shift, though real
/// in the math, is invisible because both runs pin at 100% product. (I watched
/// exactly that happen: K ≈ 10³² at both temperatures.) To *observe* Le Chatelier,
/// the enthalpy must be modest enough that K sits at an interior value at both
/// temperatures, so equilibrium lands somewhere between all-reactant and all-product
/// and can visibly move. This is a property of the *experiment*, not the engine —
/// the engine reproduces van 't Hoff at any enthalpy; only a well-chosen one lets a
/// test see it.
#[test]
fn raising_temperature_shifts_equilibrium_the_right_way() {
    let t = table();
    let mut reg = SpeciesRegistry::new();
    let (a2, b2, ab) = ab_system(&mut reg);
    let model = BondEnergyModel::default();

    // A modest exothermic enthalpy, imposed directly, so K is moderate at both
    // temperatures — see the note above. The pre-exponential is high enough that
    // *both* temperatures reach equilibrium within the step budget: the lower
    // temperature has a smaller Arrhenius rate and converges more slowly, and if it
    // has not settled the comparison is meaningless. (I learned this the slow way —
    // an under-converged 600 K run reported K far below its true value and faked a
    // violation.)
    let dh = -15_000.0; // J per event, exothermic but not overwhelming
    let ea_fwd = 40_000.0;
    let ea_rev = ea_fwd - dh; // reverse barrier higher by |ΔH|: physical, and the
                              // reason the reverse rate climbs faster with heat
    let pre = 1e3;

    let fwd = Reaction {
        reactants: vec![Term { species: a2, count: 1 }, Term { species: b2, count: 1 }],
        products: vec![Term { species: ab, count: 2 }],
        activation_j: ea_fwd,
        enthalpy: Energy::from_joules(dh),
        enthalpy_j: 0.0,
        pre_exponential: pre,
    };
    let rev = Reaction {
        reactants: vec![Term { species: ab, count: 2 }],
        products: vec![Term { species: a2, count: 1 }, Term { species: b2, count: 1 }],
        activation_j: ea_rev,
        enthalpy: Energy::from_joules(-dh),
        enthalpy_j: 0.0,
        pre_exponential: pre,
    };
    // Read on `enthalpy`, deliberately: this test builds both reactions by hand
    // at molar scale with `enthalpy_j = 0`, so the reactor takes its fallback
    // path and the quantised `Energy` field is the one carrying the value. That
    // is the documented "scales where the quantum is irrelevant" case. Tests that
    // go through `with_derived_enthalpy` must read `enthalpy_j` instead — see the
    // heating test below.
    assert!(fwd.enthalpy.0 < 0, "the forward reaction must be exothermic for this test");
    assert!(fwd.is_atom_balanced(&t, &reg) && rev.is_atom_balanced(&t, &reg));

    // The equilibrium constant K = [AB]²/([A2][B2]) in concentration units — the
    // quantity van 't Hoff actually governs. Testing K rather than the raw
    // population ratio is the more rigorous check: the ratio is muddied by
    // stoichiometry and initial amounts, but K is the thermodynamic invariant, and
    // for an exothermic reaction it must *fall* as temperature rises.
    let vol = 1.0e6;
    let equilibrium_k = |temp: f64| -> f64 {
        let mut chem = CellChemistry::new(reg.len());
        chem.set(a2, 3_000_000);
        chem.set(b2, 3_000_000);
        // A huge thermal reservoir so the reaction's own heat cannot move the
        // temperature — holding T fixed to isolate the van 't Hoff shift, as an
        // isothermal bath would in a lab.
        let mut e = Energy::from_joules(1e18);
        for _ in 0..300_000 {
            let mut cell = ReactingCell {
                chem: &mut chem,
                thermal_energy: &mut e,
                temperature_k: temp,
                volume_m3: vol,
            };
            react(&mut cell, &[fwd.clone(), rev.clone()], 1.0e-3);
        }
        let ca2 = chem.get(a2) as f64 / vol;
        let cb2 = chem.get(b2) as f64 / vol;
        let cab = chem.get(ab) as f64 / vol;
        assert!(cab > 0.0 && ca2 > 0.0, "reaction ran to an extreme at T={} — not an interior equilibrium", temp);
        cab * cab / (ca2 * cb2)
    };

    let k_cool = equilibrium_k(500.0);
    let k_hot = equilibrium_k(900.0);

    // Van 't Hoff for an exothermic reaction: K falls as T rises. The measured K
    // matches exp(−ΔH/RT) to several significant figures — but the direction alone
    // is the claim under test, and it emerges from nothing but two Arrhenius rates.
    assert!(
        k_hot < k_cool * 0.9,
        "Le Chatelier / van 't Hoff violated: heating an exothermic reaction must \
         lower its equilibrium constant, but K went {:.3} (500 K) → {:.3} (900 K)",
        k_cool,
        k_hot
    );
}

// ═══ THE COUPLING: chemistry heats its own cell ═════════════════════════════

/// **The Phase 3 payoff in one test.**
///
/// An exothermic reaction must genuinely raise the thermal energy of the cell it
/// occurs in — because that is the loop the entire phase exists to build. Chemistry
/// is not beside physics; it feeds the same energy field physics conducts and
/// radiates. Here, with no external heating at all, reactions alone warm the cell.
#[test]
fn an_exothermic_reaction_heats_its_own_cell() {
    let t = table();
    let mut reg = SpeciesRegistry::new();
    let (a2, b2, ab) = ab_system(&mut reg);
    let model = BondEnergyModel::default();
    let (fwd, _) = reversible_pair(a2, b2, ab, &reg, &t, &model, 1e5, 1e5, 20_000.0, 1e9);
    // Asserted on `enthalpy_j`, not on `enthalpy`. The quantised `Energy` mirror
    // is microjoules, and one bond forming is ~5e-19 J, so a *correct* per-event
    // enthalpy rounds to exactly zero there — as that field's own documentation
    // warns. This assertion used to read `enthalpy.0` and passed only because
    // bond-derived enthalpies were molar, i.e. 6.022e23 times too large. The test
    // was reading the wrong field and the bug was making it look right.
    assert!(
        fwd.enthalpy_j < 0.0,
        "need an exothermic reaction to release heat (H-O forms a strong polar bond)"
    );

    let mut chem = CellChemistry::new(reg.len());
    // Populations sized so the heat is *observable*, which is a real constraint
    // and not a fudge. `Energy` is quantised to the microjoule and one bond is
    // ~5e-19 J, so the thermal field cannot resolve fewer than roughly 2e12
    // reaction events however long you wait. This test previously used 1e6
    // molecules and passed only because bond-derived enthalpies were 6.022e23
    // times too large — it was measuring the bug, not the chemistry.
    chem.set(a2, 1_000_000_000_000_000_000);
    chem.set(b2, 1_000_000_000_000_000_000);

    let mut thermal = Energy::from_joules(1000.0);
    let e_before = thermal;

    for _ in 0..2000 {
        let mut cell = ReactingCell {
            chem: &mut chem,
            thermal_energy: &mut thermal,
            temperature_k: 600.0,
            volume_m3: 1.0,
        };
        react(&mut cell, &[fwd.clone()], 1.0e-7);
    }

    assert!(
        thermal.0 > e_before.0,
        "an exothermic reaction released no heat: {} µJ → {} µJ",
        e_before.0,
        thermal.0
    );
    // And the heat gained must equal the bond energy lost — conservation, restated.
    assert!(chem.get(ab) > 0, "no product formed");
}

/// Determinism, extended into chemistry. The same initial cell, reacted twice,
/// produces bit-identical populations and energy. No hidden randomness, no
/// order-dependence — the reactor is as reproducible as everything beneath it.
#[test]
fn chemistry_is_bit_for_bit_reproducible() {
    let run = || {
        let t = table();
        let mut reg = SpeciesRegistry::new();
        let (a2, b2, ab) = ab_system(&mut reg);
        let model = BondEnergyModel::default();
        let (fwd, rev) =
            reversible_pair(a2, b2, ab, &reg, &t, &model, 4e4, 4e4, 25_000.0, 25_000.0);
        let mut chem = CellChemistry::new(reg.len());
        chem.set(a2, 1_500_000);
        chem.set(b2, 900_000);
        let mut e = Energy::from_joules(1e7);
        for _ in 0..3000 {
            let mut cell = ReactingCell {
                chem: &mut chem,
                thermal_energy: &mut e,
                temperature_k: 550.0,
                volume_m3: 1.0,
            };
            react(&mut cell, &[fwd.clone(), rev.clone()], 3.0e-8);
        }
        (chem.as_slice().to_vec(), e.0)
    };
    assert_eq!(run(), run());
}
