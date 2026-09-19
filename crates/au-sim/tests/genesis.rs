//! Genesis — an open chemistry running inside the world. Phase 5a's integration.
//!
//! au-chem's network tests proved discovery in a test tube: the closure of H and
//! O atoms contains water, and generated reactions obey Boltzmann. These tests
//! prove it *in the world*: a simulation whose config declares nothing but two
//! kinds of atom, in which molecules come to exist because reactions ran in the
//! scheduler — and in which discovery is deterministic, is logged as history,
//! and survives save/resume by name.
//!
//! Two things are deliberately modest here. The heat of these reactions is
//! physically honest — one molecular event releases ~10⁻¹⁹ J — so the shared
//! thermal field barely moves; Phase 3b already proved the heat coupling at
//! scale, and genesis is about *matter*, not warmth. And the world is one cell:
//! transport of discovered species between cells is Phase 5b's problem, stated
//! now so nobody mistakes its absence for an accident.

use au_core::hash::WorldHash;
use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, Lod};
use au_physics::{energy_for_temperature, Energy, GridSpec, Mass};
use au_sim::chem_columns::{install as install_chem, species_column};
use au_sim::physics_columns::{install as install_phys, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};
use au_sim::{Config, Simulation, DEFAULT_CONFIG};

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");

const POOL: usize = 32;

/// One hot cell holding bare atoms — hydrogen at species 0, oxygen at species 1,
/// which is the order the seeds are declared in config. Everything else must be
/// reached, not placed.
struct AtomSoupGen {
    spec: GridSpec,
    pool: usize,
    temp_k: f64,
    n_h: i128,
    n_o: i128,
    materials: au_physics::MaterialRegistry,
}

impl ChunkGenerator for AtomSoupGen {
    fn generate(&self, coord: ChunkCoord, _rng: &mut au_core::Rng) -> Chunk {
        let mut c = Chunk::new(coord, Lod::Full);
        let n = self.spec.cells();
        install_phys(&mut c.columns, n);
        install_chem(&mut c.columns, n, self.pool);
        let id = self.materials.by_name("silicate").unwrap();
        let mat = self.materials.get(id).unwrap();
        let kg = 1000.0 * self.spec.cell_volume();
        let e = energy_for_temperature(kg, self.temp_k, mat, 0.0);
        c.columns.get_mut::<i128>(PHYS_MASS).unwrap().as_mut_slice().fill(Mass::from_kg(kg).0);
        c.columns
            .get_mut::<i128>(PHYS_ENERGY)
            .unwrap()
            .as_mut_slice()
            .fill(Energy::from_joules(e).0);
        c.columns.get_mut::<u16>(PHYS_MATERIAL).unwrap().as_mut_slice().fill(id.0);
        c.columns.get_mut::<i128>(species_column(0)).unwrap().as_mut_slice().fill(self.n_h);
        c.columns.get_mut::<i128>(species_column(1)).unwrap().as_mut_slice().fill(self.n_o);
        c
    }
}

fn genesis_config(seed: u64, pool: usize) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("world.seed", &seed.to_string());
    c.set("clock.seconds_per_tick", "1.0");
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", "1");
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", "1.0");
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("physics.gravity_m_s2", "0.0");

    // The whole declared chemistry: two kinds of atom. No molecules. No
    // reactions. Everything past this line has to be *found*.
    c.set("chem.enabled", "true");
    c.set("chem.open_ended", "true");
    c.set("chem.max_species", &pool.to_string());
    c.set("chem.max_atoms", "3");
    c.set("chem.barrier_j_mol", "50000");
    // A cubic centimetre holding ~10²¹ atoms — about 1 mol/L. With collision
    // prefactors derived from physics rather than fitted, concentration is no
    // longer a free choice: at the old 10⁶ molecules per m³ this chemistry
    // cannot proceed at any barrier, because the molecules essentially never
    // meet.
    c.set("chem.cell_volume_m3", "1e-6");
    c.set("chem.element.count", "2");
    c.set("chem.element.0.z", "1");
    c.set("chem.element.0.symbol", "H");
    c.set("chem.element.0.mass_mda", "1008");
    c.set("chem.element.0.valence", "1");
    c.set("chem.element.0.electronegativity_c", "220");
    c.set("chem.element.1.z", "8");
    c.set("chem.element.1.symbol", "O");
    c.set("chem.element.1.mass_mda", "15999");
    c.set("chem.element.1.valence", "2");
    c.set("chem.element.1.electronegativity_c", "344");
    c.set("chem.species.count", "2");
    c.set("chem.species.0.atoms", "1"); // H·
    c.set("chem.species.1.atoms", "8"); // O·
    c
}

fn genesis_sim(seed: u64, pool: usize) -> Simulation {
    let cfg = genesis_config(seed, pool);
    let mut sim = Simulation::new(cfg);
    let spec = GridSpec::new(1, 1, 1, 1.0);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(AtomSoupGen {
        spec,
        pool,
        temp_k: 1500.0,
        n_h: 400_000_000_000_000_000_000,
        n_o: 200_000_000_000_000_000_000,
        materials,
    }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));
    sim
}

fn water() -> au_chem::Molecule {
    au_chem::Molecule::new(
        vec![
            au_chem::Atom { z: au_chem::Z(1) },
            au_chem::Atom { z: au_chem::Z(8) },
            au_chem::Atom { z: au_chem::Z(1) },
        ],
        vec![
            au_chem::Bond::new(0, 1, au_chem::BondOrder::Single),
            au_chem::Bond::new(1, 2, au_chem::BondOrder::Single),
        ],
    )
}

fn population(sim: &Simulation, sp: usize) -> i128 {
    sim.world
        .chunks
        .iter()
        .filter_map(|(_, c)| c.columns().get::<i128>(species_column(sp)))
        .flat_map(|col| col.as_slice().iter().copied())
        .sum()
}

/// **The Phase 5a headline, in the world.** A simulation whose config names two
/// kinds of atom and nothing else runs for a while — and afterwards, water is in
/// its registry *and in its cells*. Discovered by enumeration, then actually
/// made by reactions firing in the scheduler, its existence logged as events.
#[test]
fn the_world_discovers_water_and_makes_some() {
    let mut sim = genesis_sim(101, POOL);
    let species_before = sim.world.chem_registry.len();
    assert_eq!(species_before, 2, "the config declares exactly two seed atoms");
    let events_before = sim.world.events.total();

    sim.run(60);

    assert!(
        sim.world.chem_registry.len() > species_before,
        "the world discovered nothing"
    );
    assert!(
        sim.world.events.total() > events_before,
        "discovery should have been logged as history"
    );

    let w = sim
        .world
        .chem_registry
        .lookup(&water())
        .expect("water should be in the registry — reachable, therefore found");
    let made = population(&sim, w.0 as usize);
    assert!(
        made > 0,
        "water was discovered but never made — the reactions did not fire in the loop"
    );

    // Atoms, still exact, across every discovered species: count H and O through
    // the whole registry-weighted population.
    let (mut h, mut o) = (0i128, 0i128);
    for (id, m) in sim.world.chem_registry.iter() {
        if (id.0 as usize) < POOL {
            let pop = population(&sim, id.0 as usize);
            let f = m.formula(&/* table order: H then O */ {
                let mut t = au_chem::PeriodicTable::new();
                t.add(au_chem::Element {
                    z: au_chem::Z(1),
                    symbol: "H",
                    mass_mda: 1008,
                    valence: 1,
                    electronegativity_c: 220,
                });
                t.add(au_chem::Element {
                    z: au_chem::Z(8),
                    symbol: "O",
                    mass_mda: 15999,
                    valence: 2,
                    electronegativity_c: 344,
                });
                t
            });
            h += f[0] * pop;
            o += f[1] * pop;
        }
    }
    assert_eq!(h, 400_000_000_000_000_000_000, "hydrogen atoms leaked during genesis");
    assert_eq!(o, 200_000_000_000_000_000_000, "oxygen atoms leaked during genesis");
}

/// Two worlds from the same seed discover the same molecules in the same order
/// and end bit-identical. Discovery is a pure function of history — the world
/// hash, which now folds in the registry, says so.
#[test]
fn genesis_is_deterministic() {
    let run = || {
        let mut sim = genesis_sim(202, POOL);
        sim.run(50);
        (sim.world.chem_registry.len(), sim.world.world_hash())
    };
    assert_eq!(run(), run(), "the same seed produced two different geneses");
}

/// A world saved mid-discovery, reloaded, and continued equals one that never
/// stopped — and the resumed world can *name* every species its cells hold,
/// because the registry rides in the save frame rather than being re-derived.
#[test]
fn genesis_survives_save_and_resume() {
    let mut sim = genesis_sim(303, POOL);
    sim.run(30);
    let named_at_save = sim.world.chem_registry.len();
    assert!(named_at_save > 2, "discovery should have begun before the save");
    let bytes = sim.save();

    sim.run(30);
    let straight = sim.world.world_hash();

    let mut resumed = Simulation::load(&bytes, genesis_config(303, POOL)).expect("resume");
    assert_eq!(
        resumed.world.chem_registry.len(),
        named_at_save,
        "the resumed world lost the names of its discoveries"
    );
    let spec = GridSpec::new(1, 1, 1, 1.0);
    let materials = resumed.world.materials.clone();
    resumed.world.set_generator(Box::new(AtomSoupGen {
        spec,
        pool: POOL,
        temp_k: 1500.0,
        n_h: 400_000_000_000_000_000_000,
        n_o: 200_000_000_000_000_000_000,
        materials,
    }));
    resumed.run(30);

    assert_eq!(
        straight,
        resumed.world.world_hash(),
        "a world resumed mid-discovery diverged from one that never stopped"
    );
}

/// The pool is a wall, not a crash. With shelf space for only one discovery, the
/// world keeps running, keeps its determinism, and simply never stocks what it
/// cannot shelve — populations exist only inside the pool, by construction.
#[test]
fn a_full_pool_is_a_wall_not_a_crash() {
    let run = || {
        let mut sim = genesis_sim(404, 3); // 2 seeds + room for exactly one molecule
        sim.run(40);
        sim.world.world_hash()
    };
    // Surviving the run at all is half the test; identical hashes are the rest.
    assert_eq!(run(), run(), "a cramped world lost its determinism");
}
