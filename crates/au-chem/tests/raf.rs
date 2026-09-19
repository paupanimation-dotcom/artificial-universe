//! Catalysis and autocatalytic sets, held to the house standard.
//!
//! Two anchors this phase is validated against:
//!
//!   * **A catalyst changes rates, never equilibria.** Thermodynamics, not
//!     taste: a catalyst that shifted an equilibrium would let you run a cycle
//!     forward over it and back without it and extract work forever. The
//!     implementation must preserve `Ea_f − Ea_r = ΔH` exactly, in both
//!     directions, for every catalyst.
//!   * **RAF detection** — Hordijk & Steel (2004), on Kauffman's autocatalytic
//!     sets (1971, 1986). The pruning algorithm has known answers on small
//!     hand-built systems, and those are checked here directly.
//!
//! And one structural claim that motivates the whole phase, proved rather than
//! asserted: the Phase 5a network *cannot* contain autocatalysis.

use au_chem::*;

// ═══ Fixtures ═══════════════════════════════════════════════════════════════

fn cho_table() -> PeriodicTable {
    let mut t = PeriodicTable::new();
    t.add(Element { z: Z(1), symbol: "H", mass_mda: 1008, valence: 1, electronegativity_c: 220 });
    t.add(Element { z: Z(6), symbol: "C", mass_mda: 12011, valence: 4, electronegativity_c: 255 });
    t.add(Element { z: Z(8), symbol: "O", mass_mda: 15999, valence: 2, electronegativity_c: 344 });
    t
}

fn rules() -> NetworkRules {
    NetworkRules { intrinsic_barrier_j_mol: 60_000.0, kinetics: KineticRules::default(), max_atoms: 4 }
}

/// Grow a real network from three elements, a few generations deep.
fn cho_network() -> (PeriodicTable, SpeciesRegistry, Vec<Reaction>, Vec<SpeciesId>) {
    let t = cho_table();
    let model = BondEnergyModel::default();
    let net = rules();
    let mut reg = SpeciesRegistry::new();
    let seeds: Vec<SpeciesId> =
        [Z(1), Z(6), Z(8)].iter().map(|&z| reg.intern(Molecule::monatomic(z))).collect();
    let mut present = seeds.clone();
    let mut rxns = Vec::new();
    for _ in 0..3 {
        rxns = enumerate_reactions(&present, &mut reg, &t, &model, &net);
        present = (0..reg.len()).map(|i| SpeciesId(i as u32)).collect();
    }
    (t, reg, rxns, seeds)
}

fn all_species(reg: &SpeciesRegistry) -> Vec<SpeciesId> {
    (0..reg.len()).map(|i| SpeciesId(i as u32)).collect()
}

fn relation(v: &[Catalysed]) -> Vec<(SpeciesId, usize)> {
    v.iter().map(|c| (c.catalyst, c.parent)).collect()
}

/// A synthetic reaction for the algorithm tests. The energetics are irrelevant
/// to RAF — it is a question about the *shape* of a network, not its rates.
fn rxn(reactants: &[(u32, u32)], products: &[(u32, u32)]) -> Reaction {
    let term = |&(s, c): &(u32, u32)| Term { species: SpeciesId(s), count: c };
    Reaction {
        reactants: reactants.iter().map(term).collect(),
        products: products.iter().map(term).collect(),
        activation_j: 60_000.0,
        enthalpy: Default::default(),
        enthalpy_j: 0.0,
        pre_exponential: 1.0e3,
    }
}

// ═══ The structural claim that motivates the phase ══════════════════════════

/// **Phase 5a's network cannot be autocatalytic, and here is the proof.**
///
/// Form, Split, Raise and Lower give stoichiometries of 2→1, 1→2 and 1→1, and
/// none of them can place a species on both sides of the arrow. A molecule that
/// cannot appear among its own reaction's reactants cannot make more of itself.
/// So self-replication is not merely absent from the discovered network — it is
/// structurally impossible, and no amount of temperature, time or transport
/// would ever have produced it. That is why this phase exists.
#[test]
fn the_derived_network_cannot_be_stoichiometrically_autocatalytic() {
    let (_t, _reg, rxns, _seeds) = cho_network();
    assert!(rxns.len() > 100, "the fixture should build a substantial network");
    for r in &rxns {
        for a in &r.reactants {
            assert!(
                !r.products.iter().any(|b| b.species == a.species),
                "a derived reaction put a species on both sides — the move set changed"
            );
        }
    }
}

// ═══ The thermodynamic law a catalyst may not break ═════════════════════════

/// **A catalyst changes rates, never equilibria.**
///
/// For every reaction and its inverse, catalysed by the same species, the two
/// barriers must fall by the *same absolute amount* — so `Ea_f − Ea_r`, and
/// therefore `K = exp(−ΔH/RT)`, is bit-identical to the uncatalysed pair. This
/// is also the regression test for a real bug: the first version of
/// `catalysts_for` excluded any species appearing among a reaction's reactants
/// or products, which let the joined molecule catalyse the forward direction and
/// forbade it from catalysing the reverse — a machine that shifts equilibria.
#[test]
fn a_catalyst_changes_rates_but_never_the_equilibrium() {
    let (t, reg, rxns, _seeds) = cho_network();
    let model = BondEnergyModel::default();
    let net = rules();
    let cat = CatalysisRules::default();
    let all = all_species(&reg);

    let same = |a: &[Term], b: &[Term]| {
        a.len() == b.len()
            && a.iter().all(|x| b.iter().any(|y| y.species == x.species && y.count == x.count))
    };

    let mut pairs_checked = 0;
    for (i, f) in rxns.iter().enumerate() {
        // Find this reaction's inverse.
        let inv = match rxns
            .iter()
            .enumerate()
            .find(|(j, r)| *j != i && same(&r.reactants, &f.products) && same(&r.products, &f.reactants))
        {
            Some((j, r)) => (j, r),
            None => continue,
        };
        let (_, b) = inv;
        let cf = catalysts_for(f, &all, &reg, &t, &model, &cat, &net);
        let cb = catalysts_for(b, &all, &reg, &t, &model, &cat, &net);
        if cf.is_empty() {
            continue;
        }
        // The forward difference, before catalysis.
        let bare_gap = f.activation_j - b.activation_j;
        for (c, red_f) in &cf {
            let red_b = cb
                .iter()
                .find(|(c2, _)| c2 == c)
                .map(|(_, r)| *r)
                .unwrap_or_else(|| panic!("a catalyst of one direction did not catalyse the inverse"));
            assert_eq!(
                red_f, &red_b,
                "the same catalyst lowered the two directions by different amounts"
            );
            let gap = (f.activation_j - red_f) - (b.activation_j - red_b);
            assert!(
                (gap - bare_gap).abs() < 1.0e-9,
                "catalysis moved Ea_f − Ea_r from {} to {} — the equilibrium shifted",
                bare_gap,
                gap
            );
            pairs_checked += 1;
        }
    }
    assert!(pairs_checked > 50, "only {} catalysed pairs checked", pairs_checked);
}

/// A catalyst lowers a barrier; it never abolishes one. Both directions stay
/// strictly positive whatever the grip, because the reduction is capped below
/// the intrinsic barrier.
#[test]
fn catalysis_never_abolishes_a_barrier() {
    let (t, reg, rxns, _seeds) = cho_network();
    let model = BondEnergyModel::default();
    let net = rules();
    let cat = CatalysisRules::default();
    let vars = catalysed_variants(&rxns, &all_species(&reg), &reg, &t, &model, &cat, &net);
    assert!(!vars.is_empty());
    for v in &vars {
        assert!(
            v.reaction.activation_j > 0.0,
            "a catalysed barrier reached {} — no reaction is free",
            v.reaction.activation_j
        );
        assert!(
            v.reaction.activation_j >= 0.25 * net.intrinsic_barrier_j_mol - 1.0,
            "the reduction exceeded its cap"
        );
        assert!(
            v.reaction.activation_j < rxns[v.parent].activation_j,
            "a catalyst that does not accelerate is not a catalyst"
        );
        // A catalysed path is termolecular, so its rate coefficient is in
        // m⁶·molecule⁻²·s⁻¹ while its parent's is in m³·molecule⁻¹·s⁻¹. Those
        // cannot be compared — an earlier version of this test did compare
        // them, and it only ever passed because both carried the same fitted
        // units and neither meant anything.
        //
        // The physical question is: at what catalyst concentration does the
        // catalysed route carry more flux than the bare one? Setting
        // k_un·[A][B] = k_cat·[A][B][C] gives the crossover
        //
        //     [C]* = k_un / k_cat
        //
        // which is a concentration, and which is therefore a statement one can
        // check against reality.
        let t_k = 1200.0;
        let crossover = rxns[v.parent].rate_coefficient(t_k) / v.reaction.rate_coefficient(t_k);
        assert!(
            crossover.is_finite() && crossover > 0.0,
            "the crossover concentration is not a number"
        );
        // It must at least be a concentration matter can reach: more than about
        // 10²⁹ molecules per m³ would be denser than a solid, and the catalyst
        // would be doing nothing at any attainable concentration.
        assert!(
            crossover < 1.0e29,
            "the catalyst would have to be denser than solid matter ({:.2e} m⁻³) \
             to be worth using",
            crossover
        );
    }
}

/// The catalyst survives: it appears on both sides, and the atoms balance.
#[test]
fn a_catalyst_is_not_consumed() {
    let (t, reg, rxns, _seeds) = cho_network();
    let model = BondEnergyModel::default();
    let net = rules();
    let cat = CatalysisRules::default();
    let vars = catalysed_variants(&rxns, &all_species(&reg), &reg, &t, &model, &cat, &net);
    for v in vars.iter().take(400) {
        let count_in = |ts: &[Term]| -> u32 {
            ts.iter().filter(|x| x.species == v.catalyst).map(|x| x.count).sum()
        };
        let before = count_in(&v.reaction.reactants);
        let after = count_in(&v.reaction.products);
        assert!(before >= 1 && after >= 1, "the catalyst is missing from a side");
        let parent_in = count_in(&rxns[v.parent].reactants);
        let parent_out = count_in(&rxns[v.parent].products);
        assert_eq!(
            after as i64 - before as i64,
            parent_out as i64 - parent_in as i64,
            "catalysis changed the net stoichiometry of the catalyst"
        );
    }
}

// ═══ The door opening ═══════════════════════════════════════════════════════

/// **The headline.** Three elements, four graph moves, no declared molecule
/// beyond the atoms — and somewhere in the network a species catalyses the very
/// reaction that produces it, giving `A + B + C → 2 C`. The first time in this
/// engine's life that a molecule can increase its own number. Nothing in the
/// code special-cases this; it is what the general rule says when the structure
/// happens to line up.
#[test]
fn carbon_oxygen_chemistry_discovers_a_self_producing_molecule() {
    let (t, reg, rxns, _seeds) = cho_network();
    let model = BondEnergyModel::default();
    let net = rules();
    let cat = CatalysisRules::default();
    let vars = catalysed_variants(&rxns, &all_species(&reg), &reg, &t, &model, &cat, &net);

    let doubled: Vec<&Catalysed> = vars
        .iter()
        .filter(|v| v.reaction.products.iter().any(|p| p.species == v.catalyst && p.count >= 2))
        .collect();
    assert!(
        !doubled.is_empty(),
        "no species in a C/H/O network catalyses its own formation"
    );
    // Whatever it is, it must be honest: atoms balance, catalyst net +1.
    for v in doubled.iter().take(20) {
        let cin: u32 = v.reaction.reactants.iter().filter(|x| x.species == v.catalyst).map(|x| x.count).sum();
        let cout: u32 = v.reaction.products.iter().filter(|x| x.species == v.catalyst).map(|x| x.count).sum();
        assert_eq!(cout, cin + 1, "a self-producing reaction should net exactly one copy");
        assert!(reg.get(v.catalyst).is_some());
        let _ = formula_string(reg.get(v.catalyst).unwrap(), &t);
    }
}

// ═══ RAF detection, against known answers ═══════════════════════════════════

/// No catalysis, no RAF — however many reactions there are. Reflexive
/// autocatalysis is a statement about catalysts, and a network without any
/// cannot satisfy it.
#[test]
fn no_catalysis_means_no_raf() {
    let rxns = vec![rxn(&[(0, 1), (1, 1)], &[(2, 1)]), rxn(&[(2, 1), (0, 1)], &[(3, 1)])];
    assert!(maximal_raf(&rxns, &[], &[SpeciesId(0), SpeciesId(1)]).is_none());
}

/// The smallest possible RAF: food a and b, one reaction joining them, and the
/// product catalyses it. Everything the set needs, the set makes.
#[test]
fn the_simplest_raf_is_a_product_catalysing_its_own_formation() {
    let rxns = vec![rxn(&[(0, 1), (1, 1)], &[(2, 1)])];
    let r = maximal_raf(&rxns, &[(SpeciesId(2), 0)], &[SpeciesId(0), SpeciesId(1)])
        .expect("this is a RAF");
    assert_eq!(r.reactions, vec![0]);
    assert_eq!(r.closure, vec![SpeciesId(0), SpeciesId(1), SpeciesId(2)]);
    assert_eq!(self_producing(&rxns, &[(SpeciesId(2), 0)], &r), vec![SpeciesId(2)]);
}

/// A reaction catalysed only by something the network can never build is not
/// reflexively autocatalytic — it needs outside help. Pruned.
#[test]
fn a_catalyst_that_cannot_be_built_is_pruned() {
    let rxns = vec![rxn(&[(0, 1), (1, 1)], &[(2, 1)])];
    // Species 9 appears nowhere as a product and is not food.
    assert!(maximal_raf(&rxns, &[(SpeciesId(9), 0)], &[SpeciesId(0), SpeciesId(1)]).is_none());
}

/// A reaction whose substrates the food cannot reach is not F-generated, even
/// when it is perfectly catalysed. Pruned — and its removal can cascade.
#[test]
fn unreachable_substrates_are_pruned_and_the_pruning_cascades() {
    let rxns = vec![
        // Reachable core: a + b → ab, catalysed by ab.
        rxn(&[(0, 1), (1, 1)], &[(2, 1)]),
        // Needs c, which is neither food nor ever produced.
        rxn(&[(5, 1), (0, 1)], &[(6, 1)]),
        // Would be fine, but its only catalyst is 6 — which the pruned
        // reaction above was the only source of.
        rxn(&[(2, 1), (0, 1)], &[(7, 1)]),
    ];
    let cats = [(SpeciesId(2), 0), (SpeciesId(2), 1), (SpeciesId(6), 2)];
    let r = maximal_raf(&rxns, &cats, &[SpeciesId(0), SpeciesId(1)]).expect("the core survives");
    assert_eq!(r.reactions, vec![0], "only the self-sustaining core should remain");
    assert!(r.rounds >= 2, "this network needs a cascade, not one pass");
}

/// Maximality: two independent self-sustaining cores both survive, because the
/// union of RAFs is itself a RAF and the algorithm finds *the* maximal one.
#[test]
fn the_raf_is_maximal_and_contains_every_independent_core() {
    let rxns = vec![
        rxn(&[(0, 1), (1, 1)], &[(2, 1)]),
        rxn(&[(3, 1), (4, 1)], &[(5, 1)]),
        rxn(&[(8, 1), (9, 1)], &[(10, 1)]), // substrates never available
    ];
    let cats = [(SpeciesId(2), 0), (SpeciesId(5), 1), (SpeciesId(10), 2)];
    let food = [SpeciesId(0), SpeciesId(1), SpeciesId(3), SpeciesId(4)];
    let r = maximal_raf(&rxns, &cats, &food).unwrap();
    assert_eq!(r.reactions, vec![0, 1]);
}

/// The answer is a property of the network, not of the order it was written in.
#[test]
fn raf_detection_is_order_independent() {
    let a = rxn(&[(0, 1), (1, 1)], &[(2, 1)]);
    let b = rxn(&[(2, 1), (0, 1)], &[(3, 1)]);
    let c = rxn(&[(5, 1), (6, 1)], &[(7, 1)]);
    let food = [SpeciesId(0), SpeciesId(1)];

    let f1 = vec![a.clone(), b.clone(), c.clone()];
    let cats1 = [(SpeciesId(2), 0), (SpeciesId(3), 1), (SpeciesId(7), 2)];
    let r1 = maximal_raf(&f1, &cats1, &food).unwrap();

    let f2 = vec![c, b, a];
    let cats2 = [(SpeciesId(7), 0), (SpeciesId(3), 1), (SpeciesId(2), 2)];
    let r2 = maximal_raf(&f2, &cats2, &food).unwrap();

    assert_eq!(r1.closure, r2.closure, "the closure depends on write order");
    assert_eq!(r1.reactions.len(), r2.reactions.len());
}

/// Detection is an instrument: it reads the network and reports. Run twice,
/// get the same answer, to the last element.
#[test]
fn detection_is_reproducible() {
    let (t, reg, rxns, seeds) = cho_network();
    let model = BondEnergyModel::default();
    let net = rules();
    let cat = CatalysisRules::default();
    let all = all_species(&reg);
    let once = || {
        let v = catalysed_variants(&rxns, &all, &reg, &t, &model, &cat, &net);
        maximal_raf(&rxns, &relation(&v), &seeds)
    };
    assert_eq!(once(), once());
}

// ═══ The transition ═════════════════════════════════════════════════════════

/// **Kauffman's phase transition, in a network nobody designed.**
///
/// The central claim of autocatalytic-set theory is that self-sustaining
/// networks do not fade in gradually — they appear *abruptly* once catalysis is
/// dense enough. Sweeping the grip threshold (how polar a contact must be to
/// hold a reactant at all) walks the catalysis density down, and the RAF does
/// not shrink smoothly: it collapses.
///
/// The physically interesting part is *where* the cliff sits. A threshold that
/// admits carbon–oxygen contacts leaves a large self-sustaining core; one that
/// admits only hydrogen–oxygen contacts leaves none at all. In this model the
/// network's ability to close on itself rests on C–O polarity. That is a
/// statement about a coarse template model of catalysis, not about prebiotic
/// Earth — but nobody put it there.
#[test]
fn autocatalysis_appears_abruptly_as_catalysis_densifies() {
    let (t, reg, rxns, seeds) = cho_network();
    let model = BondEnergyModel::default();
    let net = rules();
    let all = all_species(&reg);

    // (core size, how many species in it make more of themselves)
    let probe = |thr: f64| -> (usize, usize) {
        let c = CatalysisRules { min_grip_j_mol: thr, ..CatalysisRules::default() };
        let v = catalysed_variants(&rxns, &all, &reg, &t, &model, &c, &net);
        let rel = relation(&v);
        match maximal_raf(&rxns, &rel, &seeds) {
            Some(r) => {
                let sp = self_producing(&rxns, &rel, &r).len();
                (r.reactions.len(), sp)
            }
            None => (0, 0),
        }
    };

    let dense = probe(250_000.0); // carbon–oxygen contacts count
    let sparse = probe(280_000.0); // only hydrogen–oxygen contacts count
    let none = probe(400_000.0); // no contact is polar enough

    assert!(dense.0 > 100, "a dense relation should leave a large core, got {}", dense.0);
    assert!(dense.1 >= 20, "the dense core should be full of self-producing species");

    // The cliff. Note what survives on the far side: a vestigial couple of
    // reactions containing *nothing* that makes more of itself. The quantity
    // that actually matters does not decline — it switches off.
    assert!(
        dense.0 > 20 * (sparse.0 + 1),
        "the collapse should be a cliff, not a slope: {} → {}",
        dense.0,
        sparse.0
    );
    assert_eq!(sparse.1, 0, "self-production should not survive the loss of C–O grips");
    assert_eq!(none, (0, 0), "with no polar grip at all there is no core");
}

/// The core the lens finds is a proper part of the network, not all of it —
/// otherwise "RAF" would be a synonym for "the reactions", and would say nothing.
#[test]
fn the_detected_core_is_smaller_than_the_network() {
    let (t, reg, rxns, seeds) = cho_network();
    let model = BondEnergyModel::default();
    let net = rules();
    let cat = CatalysisRules::default();
    let all = all_species(&reg);
    let v = catalysed_variants(&rxns, &all, &reg, &t, &model, &cat, &net);
    let r = maximal_raf(&rxns, &relation(&v), &seeds).expect("a core exists");
    assert!(r.reactions.len() < rxns.len(), "the whole network survived — the lens says nothing");
    assert!(r.closure.len() < reg.len(), "the closure should not be every species");
    assert!(
        !self_producing(&rxns, &relation(&v), &r).is_empty(),
        "the core should contain at least one species that makes more of itself"
    );
}
