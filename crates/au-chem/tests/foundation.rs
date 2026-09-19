//! The foundation of chemistry, tested.
//!
//! Two claims have to hold before a single reaction is written, because
//! everything above depends on them:
//!
//! 1. **Atoms are conserved exactly.** The chemical analogue of Phase 2's energy
//!    invariant, and non-negotiable for the same reason: a leak is an exploit.
//!
//! 2. **The graph representation recognises molecules correctly** — the same
//!    molecule described two different ways is *seen* as the same, and two
//!    genuinely different molecules are *seen* as different. This is what makes
//!    "no hardcoded molecules" real rather than a slogan. If it fails, the
//!    molecule census that Phases 5–7 stand on is noise.

use au_chem::*;

/// A validation fixture — real measured element properties, used to check the
/// engine against known chemistry. An instrument, not content. Phase 3's chemistry
/// is generic over which elements exist; these are simply the six that built life
/// on Earth, so that "does it reproduce real chemistry?" is a question with an
/// answer.
fn life_elements() -> PeriodicTable {
    let mut t = PeriodicTable::new();
    // z, symbol, mass×1000 (daltons), valence, electronegativity×100 (Pauling)
    t.add(Element { z: Z::H, symbol: "H", mass_mda: 1_008, valence: 1, electronegativity_c: 220 });
    t.add(Element { z: Z::C, symbol: "C", mass_mda: 12_011, valence: 4, electronegativity_c: 255 });
    t.add(Element { z: Z::N, symbol: "N", mass_mda: 14_007, valence: 3, electronegativity_c: 304 });
    t.add(Element { z: Z::O, symbol: "O", mass_mda: 15_999, valence: 2, electronegativity_c: 344 });
    t.add(Element { z: Z::P, symbol: "P", mass_mda: 30_974, valence: 5, electronegativity_c: 219 });
    t.add(Element { z: Z::S, symbol: "S", mass_mda: 32_060, valence: 2, electronegativity_c: 258 });
    t
}

// ═══ ELEMENTS: conservation is an identity ═══════════════════════════════════

/// The chemical version of the most important test in the project.
///
/// Move atoms around a hundred thousand times, between arbitrary inventories, and
/// the grand total of each element is *exactly* unchanged. Not nearly. Zero
/// discrepancy. A single atom created from nothing is a free source of a limiting
/// element, and in a chemical world that is the single most valuable thing there
/// is — so something will evolve to make it, unless it is impossible.
#[test]
fn atoms_are_conserved_through_any_rearrangement() {
    let n = 6;
    let mut a = AtomCount::from_counts(vec![100, 50, 30, 80, 10, 20]);
    let mut b = AtomCount::zeros(n);
    let mut c = AtomCount::zeros(n);

    let grand_total: Vec<i128> = (0..n).map(|i| a.get(i)).collect();

    // A long, deterministic sequence of transfers between the three inventories.
    let mut seed = 12345u64;
    let mut nxt = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        seed >> 33
    };
    for _ in 0..100_000 {
        let elem = (nxt() % n as u64) as usize;
        let which = nxt() % 6;
        // Pull a small amount along a cycle a→b→c→a, so atoms genuinely circulate.
        let amount = (nxt() % 5) as i128;
        match which {
            0 if a.get(elem) >= amount => a.transfer_to(&mut b, elem, amount),
            1 if b.get(elem) >= amount => b.transfer_to(&mut c, elem, amount),
            2 if c.get(elem) >= amount => c.transfer_to(&mut a, elem, amount),
            3 if b.get(elem) >= amount => b.transfer_to(&mut a, elem, amount),
            4 if c.get(elem) >= amount => c.transfer_to(&mut b, elem, amount),
            5 if a.get(elem) >= amount => a.transfer_to(&mut c, elem, amount),
            _ => {}
        }
    }

    for i in 0..n {
        assert_eq!(
            a.get(i) + b.get(i) + c.get(i),
            grand_total[i],
            "ATOM LEAK in element index {}: {} atoms appeared or vanished",
            i,
            (a.get(i) + b.get(i) + c.get(i)) - grand_total[i]
        );
    }
}

/// Mass balances to the milli-dalton, because it is Σ count·mass with both factors
/// integer. This is the bridge to the Phase 2 mass field: when chemistry moves
/// atoms, the mass it moves is exact, so the two conservation systems never
/// disagree.
#[test]
fn mass_is_the_exact_sum_of_atomic_masses() {
    let t = life_elements();
    // One glucose molecule's worth of atoms: C6 H12 O6.
    let mut inv = AtomCount::zeros(t.len());
    inv.set(t.index_of(Z::C).unwrap(), 6);
    inv.set(t.index_of(Z::H).unwrap(), 12);
    inv.set(t.index_of(Z::O).unwrap(), 6);

    // 6·12.011 + 12·1.008 + 6·15.999 = 180.156 Da, i.e. 180_156 milli-Da. Exact.
    assert_eq!(inv.mass_mda(&t), 180_156);
}

// ═══ MOLECULES: the graph knows what it is ═══════════════════════════════════

/// Build water the way the engine will: two hydrogens, one oxygen, two O–H bonds.
/// Nobody wrote "water" anywhere. It is a graph.
fn water(t: &PeriodicTable) -> Molecule {
    let _ = t;
    Molecule::new(
        vec![Atom { z: Z::O }, Atom { z: Z::H }, Atom { z: Z::H }],
        vec![Bond::new(0, 1, BondOrder::Single), Bond::new(0, 2, BondOrder::Single)],
    )
}

/// **The test that makes "no hardcoded molecules" real.**
///
/// The same molecule, described with its atoms in a different order, must be
/// recognised as *the same molecule*. Here is water built two ways — oxygen first,
/// then hydrogen first — and the engine must see through the relabelling. If it
/// cannot, it will believe every reaction discovers a brand-new substance, and the
/// molecule census becomes meaningless.
#[test]
fn the_same_molecule_described_differently_is_recognised_as_one() {
    // Water, oxygen listed first.
    let w1 = Molecule::new(
        vec![Atom { z: Z::O }, Atom { z: Z::H }, Atom { z: Z::H }],
        vec![Bond::new(0, 1, BondOrder::Single), Bond::new(0, 2, BondOrder::Single)],
    );
    // The same water, hydrogens listed first, oxygen last — atoms permuted.
    let w2 = Molecule::new(
        vec![Atom { z: Z::H }, Atom { z: Z::H }, Atom { z: Z::O }],
        vec![Bond::new(2, 0, BondOrder::Single), Bond::new(2, 1, BondOrder::Single)],
    );

    assert_eq!(
        w1.canonical(),
        w2.canonical(),
        "the engine sees two different substances where there is only water"
    );
}

/// And it must not over-merge: two genuinely different molecules with the *same
/// formula* must be told apart. Ethanol and dimethyl ether are both C₂H₆O — same
/// atoms, different connectivity — and they are different substances with different
/// everything. Structure, not formula, is identity.
#[test]
fn isomers_with_the_same_formula_are_distinguished() {
    // Ethanol: C–C–O with hydrogens.  CH3-CH2-OH
    // atoms: 0=C 1=C 2=O 3..5=H(on C0) 6,7=H(on C1) 8=H(on O)
    let ethanol = Molecule::new(
        vec![
            Atom { z: Z::C },
            Atom { z: Z::C },
            Atom { z: Z::O },
            Atom { z: Z::H },
            Atom { z: Z::H },
            Atom { z: Z::H },
            Atom { z: Z::H },
            Atom { z: Z::H },
            Atom { z: Z::H },
        ],
        vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Single),
            Bond::new(0, 3, BondOrder::Single),
            Bond::new(0, 4, BondOrder::Single),
            Bond::new(0, 5, BondOrder::Single),
            Bond::new(1, 6, BondOrder::Single),
            Bond::new(1, 7, BondOrder::Single),
            Bond::new(2, 8, BondOrder::Single),
        ],
    );
    // Dimethyl ether: C–O–C with hydrogens.  CH3-O-CH3
    let dme = Molecule::new(
        vec![
            Atom { z: Z::C },
            Atom { z: Z::O },
            Atom { z: Z::C },
            Atom { z: Z::H },
            Atom { z: Z::H },
            Atom { z: Z::H },
            Atom { z: Z::H },
            Atom { z: Z::H },
            Atom { z: Z::H },
        ],
        vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Single),
            Bond::new(0, 3, BondOrder::Single),
            Bond::new(0, 4, BondOrder::Single),
            Bond::new(0, 5, BondOrder::Single),
            Bond::new(2, 6, BondOrder::Single),
            Bond::new(2, 7, BondOrder::Single),
            Bond::new(2, 8, BondOrder::Single),
        ],
    );

    let t = life_elements();
    // Same formula...
    assert_eq!(ethanol.formula(&t), dme.formula(&t), "these should share a formula");
    // ...different molecule.
    assert_ne!(
        ethanol.canonical(),
        dme.canonical(),
        "ethanol and dimethyl ether are different substances; the engine merged them"
    );
}

/// Bond *order* is part of identity. Ethane (C–C) and ethene (C=C) differ only in
/// how many electron pairs the carbons share, and they are different molecules with
/// different chemistry. A representation blind to bond order would think adding a
/// double bond changes nothing.
#[test]
fn bond_order_is_part_of_a_molecules_identity() {
    let ethane = Molecule::new(
        vec![Atom { z: Z::C }, Atom { z: Z::C }],
        vec![Bond::new(0, 1, BondOrder::Single)],
    );
    let ethene = Molecule::new(
        vec![Atom { z: Z::C }, Atom { z: Z::C }],
        vec![Bond::new(0, 1, BondOrder::Double)],
    );
    assert_ne!(ethane.canonical(), ethene.canonical(), "single and double bonds must differ");
}

/// A symmetric molecule stresses the canonical labeller hardest: every atom looks
/// like every other, so the tie-breaking has to be deterministic or the "canonical"
/// form is not canonical. Benzene — six carbons in a ring, alternating bond orders
/// — described starting from two different atoms must still come out identical.
#[test]
fn a_symmetric_ring_canonicalises_consistently() {
    // Benzene as a 6-ring, Kekulé form: alternating single/double.
    let ring_from = |start: usize| {
        let atoms = vec![Atom { z: Z::C }; 6];
        let mut bonds = Vec::new();
        for k in 0..6 {
            let a = (start + k) % 6;
            let b = (start + k + 1) % 6;
            let order = if k % 2 == 0 { BondOrder::Double } else { BondOrder::Single };
            bonds.push(Bond::new(a as u16, b as u16, order));
        }
        Molecule::new(atoms, bonds)
    };
    assert_eq!(
        ring_from(0).canonical(),
        ring_from(1).canonical(),
        "the same ring described from a different starting atom canonicalised differently"
    );
}

/// Breaking a bond can split a molecule in two, and the engine must be able to tell.
/// This is the check the reaction engine will lean on: after a bond breaks, is this
/// still one molecule, or has it become two that need to be counted separately?
#[test]
fn connectivity_detects_fragmentation() {
    // H–O–H is connected.
    let t = life_elements();
    assert!(water(&t).is_connected());

    // Two separate H atoms sharing no bond: two molecules wearing one struct.
    let split = Molecule::new(vec![Atom { z: Z::H }, Atom { z: Z::H }], vec![]);
    assert!(!split.is_connected(), "two unbonded atoms are not one molecule");
}

/// Valence bookkeeping: the engine must know how many bond slots an atom is using,
/// so it can refuse to give an atom more bonds than it has. Oxygen in water uses
/// two slots (its full valence); each hydrogen uses one. This is what will stop the
/// reaction engine from building impossible over-bonded molecules.
#[test]
fn used_valence_counts_bond_slots() {
    let t = life_elements();
    let w = water(&t);
    assert_eq!(w.used_valence(0), 2, "oxygen should use both its bonds");
    assert_eq!(w.used_valence(1), 1, "each hydrogen uses one");
    // Oxygen's valence is 2 and it uses 2 — saturated, no free bonds. Correct for
    // water, and the reason water is stable and unreactive.
    assert_eq!(w.used_valence(0), t.get(Z::O).unwrap().valence);
}

/// A molecule's world-hash is its canonical form, so two isomorphic molecules hash
/// identically — the determinism guarantee reaches into chemistry. The same
/// substance forming in two cells with atoms enumerated differently must not make
/// the world hash diverge.
#[test]
fn isomorphic_molecules_hash_identically() {
    use au_core::hash::WorldHash;
    let w1 = Molecule::new(
        vec![Atom { z: Z::O }, Atom { z: Z::H }, Atom { z: Z::H }],
        vec![Bond::new(0, 1, BondOrder::Single), Bond::new(0, 2, BondOrder::Single)],
    );
    let w2 = Molecule::new(
        vec![Atom { z: Z::H }, Atom { z: Z::O }, Atom { z: Z::H }],
        vec![Bond::new(1, 0, BondOrder::Single), Bond::new(1, 2, BondOrder::Single)],
    );
    assert_eq!(w1.world_hash(), w2.world_hash());
}
