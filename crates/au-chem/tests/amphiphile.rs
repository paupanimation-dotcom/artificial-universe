//! Amphiphilicity and self-assembly, held to the house standard.
//!
//! Two anchors, both outside the engine:
//!
//!   * **Traube's rule** (1891) — each further carbon strengthens a surfactant
//!     by a constant increment. Be exact about the direction: the linear factor
//!     in tail size is *input*, in one visible place, as Graham's √m was in
//!     Phase 5b. What is validated is everything around it.
//!   * **CMC buffering** — free monomer stops rising at the critical
//!     concentration and every further molecule becomes surface. That is the
//!     property that makes a membrane something which grows.

use au_chem::*;

fn table() -> PeriodicTable {
    let mut t = PeriodicTable::new();
    t.add(Element { z: Z(1), symbol: "H", mass_mda: 1008, valence: 1, electronegativity_c: 220 });
    t.add(Element { z: Z(6), symbol: "C", mass_mda: 12011, valence: 4, electronegativity_c: 255 });
    t.add(Element { z: Z(8), symbol: "O", mass_mda: 15999, valence: 2, electronegativity_c: 344 });
    t
}
fn atom(z: u8) -> Atom { Atom { z: Z(z) } }
fn single(a: usize, b: usize) -> Bond { Bond::new(a as u16, b as u16, BondOrder::Single) }
fn degree(bonds: &[Bond], i: usize) -> usize {
    bonds.iter().filter(|b| b.a as usize == i || b.b as usize == i).count()
}

/// `CH₃(CH₂)ₙ₋₁OH` — the homologous series Traube measured.
fn alcohol(n: usize) -> Molecule {
    let mut atoms: Vec<Atom> = (0..n).map(|_| atom(6)).collect();
    let mut bonds: Vec<Bond> = (0..n.saturating_sub(1)).map(|i| single(i, i + 1)).collect();
    let o = atoms.len(); atoms.push(atom(8)); bonds.push(single(n - 1, o));
    let h = atoms.len(); atoms.push(atom(1)); bonds.push(single(o, h));
    for i in 0..n {
        for _ in degree(&bonds, i)..4 {
            let k = atoms.len(); atoms.push(atom(1)); bonds.push(single(i, k));
        }
    }
    Molecule::new(atoms, bonds)
}

fn alkane(n: usize) -> Molecule {
    let mut atoms: Vec<Atom> = (0..n).map(|_| atom(6)).collect();
    let mut bonds: Vec<Bond> = (0..n.saturating_sub(1)).map(|i| single(i, i + 1)).collect();
    for i in 0..n {
        for _ in degree(&bonds, i)..4 {
            let k = atoms.len(); atoms.push(atom(1)); bonds.push(single(i, k));
        }
    }
    Molecule::new(atoms, bonds)
}

// ═══ Nothing where there should be nothing ══════════════════════════════════

/// A molecule without polar contrast is not an amphiphile, and must score
/// exactly zero. Noise here would make every gas in the world a weak surfactant.
#[test]
fn molecules_without_polar_contrast_score_zero() {
    let t = table();
    let h2 = Molecule::new(vec![atom(1), atom(1)], vec![single(0, 1)]);
    let o2 = Molecule::new(vec![atom(8), atom(8)], vec![single(0, 1)]);
    assert_eq!(amphiphilicity(&h2, &t), 0.0);
    assert_eq!(amphiphilicity(&o2, &t), 0.0);
    assert_eq!(amphiphilicity(&alkane(2), &t), 0.0, "ethane is not a surfactant");
    assert_eq!(amphiphilicity(&Molecule::monatomic(Z(6)), &t), 0.0);
}

/// A ring has no ends: every bond lies on a cycle, so no cut separates a head
/// from a tail, and the score is zero by structure rather than special case.
#[test]
fn a_ring_has_no_ends_and_therefore_no_score() {
    let t = table();
    let ring = Molecule::new(
        vec![atom(6), atom(6), atom(8)],
        vec![single(0, 1), single(1, 2), single(0, 2)],
    );
    assert!(best_cut(&ring, &t).is_none());
    assert_eq!(amphiphilicity(&ring, &t), 0.0);
}

// ═══ The shape of a surfactant ══════════════════════════════════════════════

/// A polar head on an apolar tail scores, and the cut chosen is the chemically
/// right one: the bond joining the hydroxyl to the carbon skeleton.
#[test]
fn a_polar_head_on_an_apolar_tail_is_cut_in_the_right_place() {
    let t = table();
    let m = alcohol(3);
    let cut = best_cut(&m, &t).expect("an alcohol has bridge bonds");
    assert!(cut.score > 0.0);
    assert_eq!(cut.head_atoms, 2, "the head should be the hydroxyl");
    assert!(cut.head_density > cut.tail_density * 2.0);
    assert_eq!(cut.tail_atoms, m.atom_count() - 2);
}

/// **Traube's rule**, and the increment that calibrates everything downstream:
/// a constant 3.66 per carbon, which is why one CMC decade is 7.3 score units.
#[test]
fn traube_amphiphilicity_grows_by_a_constant_increment_per_carbon() {
    let t = table();
    let scores: Vec<f64> = (1..=6).map(|n| amphiphilicity(&alcohol(n), &t)).collect();
    for w in scores.windows(2) {
        assert!(w[1] > w[0], "a longer tail must not score lower: {:?}", scores);
    }
    let inc: Vec<f64> = scores.windows(2).map(|w| w[1] - w[0]).collect();
    for i in &inc {
        assert!(
            (i - 3.66).abs() < 0.05,
            "increment {:.3} per carbon — the CMC decade of 7.3 depends on this being 3.66",
            i
        );
    }
}

/// The head has to be a head: the same skeleton without a hydroxyl scores far
/// lower.
#[test]
fn the_same_skeleton_without_a_head_is_not_a_surfactant() {
    let t = table();
    for n in 2..=5 {
        let with_head = amphiphilicity(&alcohol(n), &t);
        let without = amphiphilicity(&alkane(n), &t);
        assert!(with_head > 10.0 * without.max(0.01), "n={}: {} vs {}", n, with_head, without);
    }
}

/// Water is polar all over rather than polar at one end — a solvent, not a
/// surfactant.
#[test]
fn water_is_a_solvent_not_a_surfactant() {
    let t = table();
    let water = Molecule::new(vec![atom(1), atom(8), atom(1)], vec![single(0, 1), single(1, 2)]);
    assert!(amphiphilicity(&alcohol(5), &t) > 20.0 * amphiphilicity(&water, &t).max(0.01));
}

/// The same molecule written with its atoms in a different order is the same
/// molecule. If the score moved it would be a property of the file format.
#[test]
fn the_score_does_not_depend_on_atom_numbering() {
    let t = table();
    let a = Molecule::new(
        vec![atom(6), atom(1), atom(1), atom(1), atom(8), atom(1)],
        vec![single(0, 1), single(0, 2), single(0, 3), single(0, 4), single(4, 5)],
    );
    let b = Molecule::new(
        vec![atom(8), atom(1), atom(6), atom(1), atom(1), atom(1)],
        vec![single(0, 1), single(0, 2), single(2, 3), single(2, 4), single(2, 5)],
    );
    assert_eq!(a.canonical(), b.canonical(), "the fixture is wrong");
    assert!((amphiphilicity(&a, &t) - amphiphilicity(&b, &t)).abs() < 1.0e-12);
}

// ═══ Self-assembly ══════════════════════════════════════════════════════════

/// **Free monomer is buffered at the CMC.** Ten times the surfactant gives the
/// same dissolved concentration and nine times the surface. This is what makes
/// membranes grow by accretion instead of becoming a stronger solution.
#[test]
fn free_monomer_is_buffered_at_the_critical_concentration() {
    let t = table();
    let r = MembraneRules::default();
    let score = amphiphilicity(&alcohol(4), &t);
    let volume = 1.0e-15; // a cubic micron
    let cap = (critical_concentration(score, &r).unwrap() * volume).floor() as i128;
    assert!(cap > 0, "a micron-scale cell must hold a measurable free population");
    for mult in [1i128, 2, 10, 100] {
        let pop = cap * mult;
        let free = pop - assembled_count(pop, score, volume, &r);
        assert_eq!(free, cap.min(pop), "free monomer left the buffer at ×{}", mult);
    }
}

/// Nothing is invented and nothing lost: free plus assembled is the population,
/// exactly, at every concentration.
#[test]
fn assembly_partitions_exactly() {
    let t = table();
    let r = MembraneRules::default();
    for n in 1..=6 {
        let score = amphiphilicity(&alcohol(n), &t);
        for pop in [0i128, 1, 999, 1_000_000_000, 999_999_999_999_999_999] {
            let a = assembled_count(pop, score, 1.0e-15, &r);
            assert!(a >= 0 && a <= pop);
        }
    }
}

/// **Traube carried through to the CMC**: a longer tail assembles at a lower
/// concentration, by decades.
#[test]
fn a_stronger_amphiphile_assembles_at_a_lower_concentration() {
    let t = table();
    let r = MembraneRules::default();
    let c1 = critical_concentration(amphiphilicity(&alcohol(1), &t), &r).unwrap();
    let c6 = critical_concentration(amphiphilicity(&alcohol(6), &t), &r).unwrap();
    assert!(c1 > c6);
    assert!(c1 / c6 > 100.0, "the CMC should fall by decades, got {:.1}", c1 / c6);
}

/// A molecule without polar contrast never assembles, however concentrated.
/// Shape makes a membrane, not crowding.
#[test]
fn a_non_amphiphile_never_assembles_at_any_concentration() {
    let r = MembraneRules::default();
    assert!(critical_concentration(0.0, &r).is_none());
    assert_eq!(assembled_count(i128::MAX / 2, 0.0, 1.0e-15, &r), 0);
}

// ═══ The geometry that sets the scale ═══════════════════════════════════════

/// **The arithmetic that says a compartment must be cell-sized.**
///
/// Closing a bag needs `enclosing_area / area_per_molecule` molecules, and with
/// a 0.4 nm² headgroup that number states the scale by itself: ~10¹⁹ for a cubic
/// metre (hopeless), ~10⁹ for a eukaryote-sized volume, ~10⁷ for a bacterium.
///
/// The last is an external check the model was never fitted to: *E. coli*'s
/// membrane holds on the order of 10⁷ lipids. Given only a headgroup area and
/// the geometry of a sphere, the arithmetic lands on the size of a real
/// bacterial membrane.
#[test]
fn closing_a_compartment_takes_a_cell_sized_number_of_molecules() {
    let r = MembraneRules::default();
    let needed = |v: f64| enclosing_area(v) / r.area_per_molecule_m2;

    let tank = needed(1.0);
    assert!((1.0e19..1.0e20).contains(&tank), "a cubic metre: {:.2e}", tank);

    let eukaryote = needed(1.0e-15);
    assert!((1.0e9..1.0e10).contains(&eukaryote), "a 10 µm cell: {:.2e}", eukaryote);

    // The anchor: a bacterium-sized compartment, against a real bacterium.
    let bacterium = needed(1.0e-18);
    assert!(
        (5.0e6..5.0e7).contains(&bacterium),
        "a 1 µm³ compartment needs {:.2e} amphiphiles to close; E. coli's membrane \
         holds on the order of 10⁷ lipids",
        bacterium
    );

    // Scaling is geometric: a thousandfold volume needs a hundredfold surface.
    assert!((eukaryote / bacterium / 100.0 - 1.0).abs() < 0.01);
}

/// Coverage is the fraction of the boundary that is actually covered, and
/// permeability follows from what is left open — never reaching zero, because a
/// real bilayer leaks and a perfect wall would make a tomb rather than a cell.
#[test]
fn coverage_sets_permeability_and_a_sealed_bag_still_leaks() {
    let r = MembraneRules::default();
    let v = 1.0e-18;
    let need = (enclosing_area(v) / r.area_per_molecule_m2) as i128;

    assert_eq!(coverage(0, v, &r), 0.0);
    assert!((coverage(need, v, &r) - 1.0).abs() < 0.01, "a full complement should just close it");
    assert!(coverage(2 * need, v, &r) > 1.9, "twice the material covers twice the boundary");

    assert!((permeability(0.0, &r) - 1.0).abs() < 1.0e-12, "an open cell is fully open");
    assert_eq!(permeability(1.0, &r), r.residual_permeability);
    assert_eq!(permeability(5.0, &r), r.residual_permeability, "coverage cannot over-seal");
    // Monotone in between.
    assert!(permeability(0.25, &r) > permeability(0.75, &r));
}

/// Applied to a vocabulary the engine invented: only some species assemble, and
/// those that do genuinely have a head and a tail.
#[test]
fn a_discovered_network_assembles_only_its_surfactants() {
    let t = table();
    let model = BondEnergyModel::default();
    let net = NetworkRules {
        intrinsic_barrier_j_mol: 60_000.0,
        kinetics: KineticRules::default(),
        max_atoms: 6,
    };
    let mut reg = SpeciesRegistry::new();
    let seeds: Vec<SpeciesId> =
        [Z(1), Z(6), Z(8)].iter().map(|&z| reg.intern(Molecule::monatomic(z))).collect();
    let mut present = seeds.clone();
    for _ in 0..3 {
        let _ = enumerate_reactions(&present, &mut reg, &t, &model, &net);
        present = (0..reg.len()).map(|i| SpeciesId(i as u32)).collect();
    }
    assert!(reg.len() > 200, "the fixture should discover a substantial vocabulary");

    let scores = score_table((0..reg.len()).map(|i| reg.get(SpeciesId(i as u32)).unwrap()), &t);
    let assemblers = scores.iter().filter(|&&s| s >= MembraneRules::default().min_score).count();
    assert!(assemblers > 0, "nothing in the network can build a membrane");
    assert!(
        assemblers < reg.len() / 2,
        "{} of {} species assembled — the edge is too generous",
        assemblers,
        reg.len()
    );
    for (i, &sc) in scores.iter().enumerate() {
        if sc >= MembraneRules::default().min_score {
            let cut = best_cut(reg.get(SpeciesId(i as u32)).unwrap(), &t).unwrap();
            assert!(cut.tail_atoms >= 2 && cut.head_density > cut.tail_density);
        }
    }
}
