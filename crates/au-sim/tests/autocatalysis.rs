//! Catalysis and autocatalysis in the running world — Phase 5c's integration.
//!
//! au-chem proved the relation and the algorithm. These tests prove the tower:
//! that catalysis is genuinely off by default (so every earlier world is
//! untouched), that catalysed reactions accelerate chemistry without breaking
//! the atom books, and that the survey is what it claims to be — an instrument
//! that changes nothing it measures.

use au_core::hash::WorldHash;
use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, Lod};
use au_physics::{energy_for_temperature, Energy, GridSpec, Mass};
use au_sim::autocatalysis::survey;
use au_sim::chem_columns::{install as install_chem, species_column};
use au_sim::physics_columns::{install as install_phys, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};
use au_sim::{Config, Simulation, DEFAULT_CONFIG};

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");
const POOL: usize = 48;

/// One hot cell holding bare atoms — the smallest world in which chemistry can
/// have a history.
struct SoupGen {
    pool: usize,
    temp_k: f64,
    seeds: Vec<(usize, i128)>,
    materials: au_physics::MaterialRegistry,
}

impl ChunkGenerator for SoupGen {
    fn generate(&self, coord: ChunkCoord, _rng: &mut au_core::Rng) -> Chunk {
        let mut c = Chunk::new(coord, Lod::Full);
        install_phys(&mut c.columns, 1);
        install_chem(&mut c.columns, 1, self.pool);
        let id = self.materials.by_name("silicate").unwrap();
        let mat = self.materials.get(id).unwrap();
        let kg = 1000.0;
        let e = energy_for_temperature(kg, self.temp_k, mat, 0.0);
        c.columns.get_mut::<i128>(PHYS_MASS).unwrap().as_mut_slice()[0] = Mass::from_kg(kg).0;
        c.columns.get_mut::<i128>(PHYS_ENERGY).unwrap().as_mut_slice()[0] =
            Energy::from_joules(e).0;
        c.columns.get_mut::<u16>(PHYS_MATERIAL).unwrap().as_mut_slice()[0] = id.0;
        for &(sp, n) in &self.seeds {
            c.columns.get_mut::<i128>(species_column(sp)).unwrap().as_mut_slice()[0] = n;
        }
        c
    }
}

/// A C/H/O world: three elements, no declared molecules, open chemistry.
fn cho_config(seed: u64, catalysis: bool) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("world.seed", &seed.to_string());
    // Nanosecond ticks, because that is how fast this chemistry actually is.
    // A ~1.7 M solution at 2500 K over a 200 kJ/mol barrier gives each molecule
    // ~2×10⁷ successful collisions per second — a lifetime of tens of
    // nanoseconds. With real collision prefactors the timescale is no longer
    // ours to choose; it is dictated by the concentration and the barrier.
    c.set("clock.seconds_per_tick", "1e-9");
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", "1");
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", "1.0");
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("physics.gravity_m_s2", "0.0");
    c.set("chem.enabled", "true");
    // One cubic centimetre, holding ~10²¹ atoms: about 1.7 mol/L, a real
    // concentrated solution.
    //
    // The previous value was 10⁶ m³ — a cube a hundred metres on a side — chosen
    // so that catalysis would be visible against a *fitted* prefactor. That
    // crutch is gone: prefactors are now derived from collision theory, so the
    // cell has to be a real cell and the concentration a real concentration.
    // Catalysis is visible here for a physical reason instead: the crossover
    // concentration at which the catalysed route overtakes is ~0.1 mol/L, and
    // this world is above it.
    c.set("chem.cell_volume_m3", "1e-6");
    c.set("chem.open_ended", "true");
    c.set("chem.max_species", &POOL.to_string());
    c.set("chem.max_atoms", "3");
    // A high intrinsic barrier. A catalysed reaction is termolecular — it must
    // find its catalyst as well as both reactants — and pays an encounter volume
    // for the privilege, so it only wins when the barrier it removes is worth
    // more than that encounter costs. Phase 5c found this empirically; the
    // kinetics module now says why, and puts a number on it.
    c.set("chem.barrier_j_mol", "200000");
    c.set("chem.element.count", "3");
    c.set("chem.element.0.z", "1");
    c.set("chem.element.0.symbol", "H");
    c.set("chem.element.0.mass_mda", "1008");
    c.set("chem.element.0.valence", "1");
    c.set("chem.element.0.electronegativity_c", "220");
    c.set("chem.element.1.z", "6");
    c.set("chem.element.1.symbol", "C");
    c.set("chem.element.1.mass_mda", "12011");
    c.set("chem.element.1.valence", "4");
    c.set("chem.element.1.electronegativity_c", "255");
    c.set("chem.element.2.z", "8");
    c.set("chem.element.2.symbol", "O");
    c.set("chem.element.2.mass_mda", "15999");
    c.set("chem.element.2.valence", "2");
    c.set("chem.element.2.electronegativity_c", "344");
    c.set("chem.species.count", "3");
    c.set("chem.species.0.atoms", "1");
    c.set("chem.species.1.atoms", "6");
    c.set("chem.species.2.atoms", "8");
    if catalysis {
        c.set("chem.catalysis", "true");
    }
    c
}

fn cho_sim(seed: u64, catalysis: bool, temp_k: f64) -> Simulation {
    let mut sim = Simulation::new(cho_config(seed, catalysis));
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(SoupGen {
        pool: POOL,
        temp_k,
        seeds: vec![(0, 600_000_000_000_000_000_000), (1, 200_000_000_000_000_000_000), (2, 300_000_000_000_000_000_000)],
        materials,
    }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));
    sim
}

/// Every species' population in the one cell.
fn pops(sim: &Simulation) -> Vec<i128> {
    let chunk = sim.world.chunks.iter().next().unwrap().1;
    (0..POOL)
        .map(|s| chunk.columns().get::<i128>(species_column(s)).unwrap().as_slice()[0])
        .collect()
}

/// Total atoms of each element across every discovered species.
fn atom_audit(sim: &Simulation) -> (i128, i128, i128) {
    let table = {
        let mut t = au_chem::PeriodicTable::new();
        t.add(au_chem::Element { z: au_chem::Z(1), symbol: "H", mass_mda: 1008, valence: 1, electronegativity_c: 220 });
        t.add(au_chem::Element { z: au_chem::Z(6), symbol: "C", mass_mda: 12011, valence: 4, electronegativity_c: 255 });
        t.add(au_chem::Element { z: au_chem::Z(8), symbol: "O", mass_mda: 15999, valence: 2, electronegativity_c: 344 });
        t
    };
    let p = pops(sim);
    let (mut h, mut c, mut o) = (0i128, 0i128, 0i128);
    for (id, m) in sim.world.chem_registry.iter() {
        let i = id.0 as usize;
        if i < POOL && p[i] != 0 {
            let f = m.formula(&table);
            h += f[0] * p[i];
            c += f[1] * p[i];
            o += f[2] * p[i];
        }
    }
    (h, c, o)
}

/// **Catalysis is off unless asked for.** An omitted flag and an explicit
/// `false` give the same world, and that world is the Phase 5b one. Turning it
/// on changes the chemistry — which is the point, and is why it must be opt-in.
#[test]
fn catalysis_is_off_by_default_and_opt_in_changes_the_world() {
    let mut a = cho_sim(31, false, 2000.0);
    a.run(40);
    let mut b = Simulation::new({
        let mut c = cho_config(31, false);
        c.set("chem.catalysis", "false");
        c
    });
    let materials = b.world.materials.clone();
    b.world.set_generator(Box::new(SoupGen {
        pool: POOL,
        temp_k: 2000.0,
        seeds: vec![(0, 600_000_000_000_000_000_000), (1, 200_000_000_000_000_000_000), (2, 300_000_000_000_000_000_000)],
        materials,
    }));
    b.world.activate(ChunkCoord::new(0, 0, 0));
    b.run(40);
    assert_eq!(a.world.world_hash(), b.world.world_hash(), "an explicit false differs from absent");

    let mut c = cho_sim(31, true, 2000.0);
    c.run(40);
    assert_ne!(
        a.world.world_hash(),
        c.world.world_hash(),
        "enabling catalysis changed nothing — the variants are not reaching the reactor"
    );
}

/// **The atom books survive catalysis.** A catalysed reaction writes the same
/// species on both sides, and the terms are coalesced (`AB + AB` becomes
/// `2 AB`). If that bookkeeping were wrong the reactor would quietly create or
/// destroy matter, so this is the test that guards it.
#[test]
fn atoms_are_conserved_with_catalysis_running() {
    let mut sim = cho_sim(32, true, 2500.0);
    let before = atom_audit(&sim);
    sim.run(300);
    let after = atom_audit(&sim);
    assert_eq!(
        before,
        (600_000_000_000_000_000_000, 200_000_000_000_000_000_000, 300_000_000_000_000_000_000),
        "the seeding is wrong"
    );
    assert_eq!(after, before, "catalysis broke the atom books");
    assert!(pops(&sim).iter().all(|&v| v >= 0), "a population went negative");
}

/// A catalyst makes chemistry *faster*. Same world, same seed, same temperature:
/// the catalysed one has travelled further from its starting atoms by the time
/// the uncatalysed one is still working on it.
#[test]
fn catalysis_accelerates_the_chemistry() {
    let mut plain = cho_sim(33, false, 1400.0);
    let mut cata = cho_sim(33, true, 1400.0);
    plain.run(60);
    cata.run(60);
    // Progress = atoms no longer sitting in the monatomic seeds.
    let free = |s: &Simulation| -> i128 { pops(s)[0] + pops(s)[1] + pops(s)[2] };
    assert!(
        free(&cata) < free(&plain),
        "catalysed {} vs plain {} free atoms — catalysis did not accelerate anything",
        free(&cata),
        free(&plain)
    );
}

/// **The instrument changes nothing it measures.** Surveying a world leaves its
/// hash, its clock and its registry exactly as they were — the enumeration runs
/// on a clone precisely so that looking cannot intern a species into the world.
#[test]
fn the_survey_does_not_disturb_the_world() {
    let mut sim = cho_sim(34, true, 2500.0);
    sim.run(120);
    let hash_before = sim.world.world_hash();
    let species_before = sim.world.chem_registry.len();

    let s = survey(&sim.world).expect("an open-chemistry world can be surveyed");
    assert!(!s.reactions.is_empty(), "the survey found no network at all");

    assert_eq!(sim.world.world_hash(), hash_before, "the survey changed the world hash");
    assert_eq!(sim.world.chem_registry.len(), species_before, "the survey interned a species");

    // And it is repeatable.
    let s2 = survey(&sim.world).unwrap();
    assert_eq!(s.raf.is_some(), s2.raf.is_some());
    assert_eq!(s.self_producing, s2.self_producing);
}

/// **A world whose chemistry has closed on itself.** Three elements, no
/// declared molecule, catalysis derived from structure — and the lens finds a
/// set of reactions that collectively makes its own catalysts from bare atoms,
/// containing at least one species that catalyses a reaction producing it.
#[test]
fn a_carbon_world_closes_on_itself() {
    let mut sim = cho_sim(35, true, 2500.0);
    sim.run(200);
    let s = survey(&sim.world).expect("surveyable");
    assert!(s.catalyst_count() > 0, "no species in this world catalyses anything");
    let raf = s.raf.as_ref().expect("this network should have an autocatalytic core");
    assert!(!raf.reactions.is_empty());
    assert!(
        raf.closure.len() >= s.food.len(),
        "the closure must at least contain the food it was given"
    );
    assert!(
        !s.self_producing.is_empty(),
        "the core contains nothing that makes more of itself"
    );
    // The core is a proper part of the network — otherwise the lens says nothing.
    assert!(raf.reactions.len() <= s.reactions.len());
}

/// Without catalysis there is no relation, and therefore no RAF — the same
/// world, looked at through the same lens, correctly reports that its chemistry
/// does not sustain itself.
#[test]
fn a_world_without_catalysis_has_no_autocatalytic_core() {
    let mut sim = cho_sim(36, false, 2500.0);
    sim.run(200);
    let s = survey(&sim.world).expect("surveyable");
    assert_eq!(s.catalyst_count(), 0, "catalysis is off; nothing should catalyse");
    assert!(s.raf.is_none(), "a network with no catalysts cannot be reflexively autocatalytic");
    assert!(!s.reactions.is_empty(), "…but there is still plenty of chemistry");
}

/// Catalysis does not cost determinism: same seed, same world, to the bit — and
/// a world saved mid-reaction resumes identically.
#[test]
fn catalysed_worlds_are_reproducible_and_resumable() {
    let run = || {
        let mut s = cho_sim(37, true, 2500.0);
        s.run(150);
        s.world.world_hash()
    };
    assert_eq!(run(), run(), "catalysis lost determinism");

    let mut sim = cho_sim(37, true, 2500.0);
    sim.run(75);
    let bytes = sim.save();
    sim.run(75);
    let straight = sim.world.world_hash();

    let mut resumed = Simulation::load(&bytes, cho_config(37, true)).expect("resume");
    let materials = resumed.world.materials.clone();
    resumed.world.set_generator(Box::new(SoupGen {
        pool: POOL,
        temp_k: 2500.0,
        seeds: vec![(0, 600_000_000_000_000_000_000), (1, 200_000_000_000_000_000_000), (2, 300_000_000_000_000_000_000)],
        materials,
    }));
    resumed.run(75);
    assert_eq!(straight, resumed.world.world_hash(), "a catalysed world diverged on resume");
}
