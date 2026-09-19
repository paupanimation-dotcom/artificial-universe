//! A planet, generating its own world — Phase 4.
//!
//! au-planet proved the parameter → state derivation (the habitable zone emerges
//! from Stefan–Boltzmann and water's phase diagram). These tests prove the
//! *generation*: a world that declares a planet in config fills its own chunks —
//! no hand-seeded cells, no test fixture generator — and the matter placed obeys
//! the physics, not the seed. The seed decides where the rock stands; the star
//! decides what everything is.

use au_core::hash::WorldHash;
use au_data::chunk::ChunkCoord;
use au_physics::{derive, Energy, GridSpec, MaterialId, AG_PER_KG};
use au_sim::physics_columns::{PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};
use au_sim::{Config, Simulation, DEFAULT_CONFIG};

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");

/// A planet-world config: reference materials, a small grid, and a declared
/// planet. `orbit_au` is the knob the key test turns — everything else identical.
fn planet_config(seed: u64, orbit_au: f64) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("world.seed", &seed.to_string());
    c.set("clock.seconds_per_tick", "1.0");
    c.set("clock.rescale.enabled", "false");

    // Chunk grid: 8 × 1 × 16 cells of 5 m — a 40 m wide, 80 m tall slice.
    c.set("physics.grid.nx", "8");
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "16");
    c.set("physics.grid.cell_m", "5.0");
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("physics.gravity_m_s2", "0.0");

    // The planet. Albedo 0 — a perfectly dark world — stands in for the
    // greenhouse warming this phase does not model. A telling detail: even at
    // albedo 0, 1 AU only reaches ~277 K; you cannot reach Earth's actual 288 K
    // surface from albedo alone, because the missing 11 K *is* the greenhouse.
    // The warm test orbit sits at 0.95 AU for a few degrees of margin; the same
    // config swings from ocean to ice purely by moving the orbit.
    c.set("planet.enabled", "true");
    c.set("planet.orbit_au", &orbit_au.to_string());
    c.set("planet.albedo", "0.0");
    c.set("planet.relief_m", "30.0");
    c.set("planet.relief_wavelength_m", "60.0");
    c
}

/// The chunk just below sea level — z spans −80 m .. 0 m — where rock, volatile,
/// and the terrain's decision between them all live.
const SHORE: ChunkCoord = ChunkCoord { x: 0, y: 0, z: -1 };

fn spec() -> GridSpec {
    GridSpec::new(8, 1, 16, 5.0)
}

/// Read one cell's (mass_ug, energy_uj, material) from a chunk in the sim.
fn cell(sim: &Simulation, coord: ChunkCoord, i: usize) -> (i128, i128, u16) {
    let ch = sim.world.chunks.get(coord).expect("chunk should exist");
    let m = ch.columns().get::<i128>(PHYS_MASS).unwrap().as_slice()[i];
    let e = ch.columns().get::<i128>(PHYS_ENERGY).unwrap().as_slice()[i];
    let mat = ch.columns().get::<u16>(PHYS_MATERIAL).unwrap().as_slice()[i];
    (m, e, mat)
}

/// **No hands.** A world that declares a planet fills its own chunks: activating a
/// coordinate produces rock and volatile cells without any test generator being
/// installed. This is the first time in the project that matter exists without a
/// fixture placing it.
#[test]
fn a_planet_world_fills_its_own_chunks() {
    let mut sim = Simulation::new(planet_config(11, 1.0));
    let coords: Vec<ChunkCoord> = (0..4).map(|x| ChunkCoord { x, y: 0, z: -1 }).collect();
    for c in &coords {
        sim.world.activate(*c);
    }
    sim.run(1);

    let s = spec();
    let mut rock = 0;
    let mut volatile = 0;
    let mut vacuum = 0;
    let mats = sim.world.materials.clone();
    let water = mats.by_name("water").unwrap();
    let silicate = mats.by_name("silicate").unwrap();
    for c in &coords {
        for i in 0..s.cells() {
            let (_, _, mat) = cell(&sim, *c, i);
            if mat == silicate.0 {
                rock += 1;
            } else if mat == water.0 {
                volatile += 1;
            } else if mat == MaterialId::VACUUM.0 {
                vacuum += 1;
            }
        }
    }
    assert!(rock > 0, "the planet placed no rock");
    assert!(volatile > 0, "below sea level, hollows should hold the volatile");
    assert_eq!(rock + volatile + vacuum, 4 * s.cells());
}

/// In every column, matter is ordered by cause: rock at the bottom (it is where
/// the terrain stands), volatile above it up to sea level, nothing floating. Going
/// up a column the sequence is rock… volatile… — never volatile under rock, never
/// a gap inside the fill. Structure follows from the fill rule, and the fill rule
/// follows from where things belong.
#[test]
fn columns_are_ordered_rock_then_volatile() {
    let mut sim = Simulation::new(planet_config(12, 1.0));
    sim.world.activate(SHORE);
    sim.run(1);

    let s = spec();
    let mats = sim.world.materials.clone();
    let water = mats.by_name("water").unwrap().0;
    let silicate = mats.by_name("silicate").unwrap().0;

    for x in 0..s.nx {
        // Rank: rock=0, volatile=1, vacuum=2 — must be non-decreasing going up.
        let mut rank = 0;
        for z in 0..s.nz {
            let (_, _, mat) = cell(&sim, SHORE, s.idx(x, 0, z));
            let r = if mat == silicate {
                0
            } else if mat == water {
                1
            } else {
                2
            };
            assert!(
                r >= rank,
                "column x={} broke order at z={}: material {} below rank {}",
                x,
                z,
                mat,
                rank
            );
            rank = r;
        }
    }
}

/// **The Phase 4 headline, in the sim.** Two worlds, identical in every declared
/// respect except orbital distance. In the near world the below-sea-level cells
/// derive as **liquid** — an ocean. In the far world the *same* cells, placed by
/// the *same* generator code, derive as **solid** — an ice sheet. The generator
/// never chose; it wrote the energy the star's temperature implies, and the phase
/// fell out of the physics when we looked.
#[test]
fn a_cold_planet_makes_ice_where_a_warm_one_makes_ocean() {
    let phase_of_volatile = |orbit_au: f64| -> au_physics::Phase {
        let mut sim = Simulation::new(planet_config(13, orbit_au));
        let coords: Vec<ChunkCoord> = (0..4).map(|x| ChunkCoord { x, y: 0, z: -1 }).collect();
        for c in &coords {
            sim.world.activate(*c);
        }
        sim.run(1);
        let s = spec();
        let mats = sim.world.materials.clone();
        let water_id = mats.by_name("water").unwrap();
        let water = mats.get(water_id).unwrap();
        for c in &coords {
            for i in 0..s.cells() {
                let (m, e, mat) = cell(&sim, *c, i);
                if mat == water_id.0 {
                    let kg = m as f64 / AG_PER_KG as f64;
                    let state = derive(kg, Energy(e).as_joules(), water, 101_325.0, s.cell_volume());
                    return state.phase;
                }
            }
        }
        panic!("no volatile cell found at {} AU", orbit_au);
    };

    assert_eq!(
        phase_of_volatile(0.95),
        au_physics::Phase::Liquid,
        "at 0.95 AU (albedo-0 greenhouse proxy) the volatile should be an ocean"
    );
    assert_eq!(
        phase_of_volatile(1.6),
        au_physics::Phase::Solid,
        "at 1.6 AU the same cells should be an ice sheet — physics decided, not the generator"
    );
}

/// The generated world is deterministic: two boots of the same seed, activating
/// the same chunks and running the same ticks, hash identically — the planet is
/// part of the world function, not a source of drift.
#[test]
fn a_planet_world_is_deterministic() {
    let run = || {
        let mut sim = Simulation::new(planet_config(14, 1.0));
        sim.world.activate(SHORE);
        sim.world.activate(ChunkCoord { x: 1, y: 0, z: -1 });
        sim.run(50);
        sim.world.world_hash()
    };
    assert_eq!(run(), run(), "the same seed generated two different planets");
}

/// Save/resume with a planet: the generator is reinstalled automatically by
/// `Simulation::load` (it is config, and config rides in the snapshot), so a
/// resumed planet world equals one that never stopped — including any chunk the
/// resumed world has to regenerate from seed.
#[test]
fn a_planet_world_survives_save_and_resume() {
    let mut sim = Simulation::new(planet_config(15, 1.0));
    sim.world.activate(SHORE);
    sim.run(40);
    let bytes = sim.save();

    sim.run(40);
    let straight = sim.world.world_hash();

    let mut resumed = Simulation::load(&bytes, planet_config(15, 1.0)).expect("resume");
    resumed.run(40);
    let rejoined = resumed.world.world_hash();

    assert_eq!(straight, rejoined, "a resumed planet diverged from one that never stopped");
}
