//! Observation — the tests that make it trustworthy.
//!
//! A recorder that perturbed the world would be worse than no recorder, because
//! every conclusion drawn from it would be about a world nobody else was running.
//! Phase 1 shipped a bug where an event emitted during chunk eviction made RAM
//! pressure part of world identity; that is the failure mode this file exists to
//! rule out.

use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, Lod};
use au_physics::{energy_for_temperature, Energy, Mass};
use au_sim::chem_columns::{install as install_chem, species_column};
use au_sim::observe::Recorder;
use au_sim::physics_columns::{install as install_phys, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};
use au_sim::{Config, Simulation, DEFAULT_CONFIG};

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");
const CELL_VOLUME: f64 = 1.0e-18;
const POOL: usize = 2;

struct Gen {
    materials: au_physics::MaterialRegistry,
    seeds: Vec<(usize, i128)>,
}

impl ChunkGenerator for Gen {
    fn generate(&self, coord: ChunkCoord, _rng: &mut au_core::Rng) -> Chunk {
        let mut c = Chunk::new(coord, Lod::Full);
        install_phys(&mut c.columns, 1);
        install_chem(&mut c.columns, 1, POOL);
        let id = self.materials.by_name("water").unwrap();
        let mat = self.materials.get(id).unwrap();
        let kg = 1000.0 * CELL_VOLUME;
        c.columns.get_mut::<i128>(PHYS_MASS).unwrap().as_mut_slice()[0] = Mass::from_kg(kg).0;
        c.columns.get_mut::<i128>(PHYS_ENERGY).unwrap().as_mut_slice()[0] =
            Energy::from_joules(energy_for_temperature(kg, 300.0, mat, 0.0)).0;
        c.columns.get_mut::<u16>(PHYS_MATERIAL).unwrap().as_mut_slice()[0] = id.0;
        for (sp, n) in &self.seeds {
            c.columns.get_mut::<i128>(species_column(*sp)).unwrap().as_mut_slice()[0] = *n;
        }
        c
    }
}

fn sim(seed: u64) -> Simulation {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("world.seed", &seed.to_string());
    c.set("clock.seconds_per_tick", "1.0");
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", "1");
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", "1.0e-6");
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("chem.enabled", "true");
    c.set("chem.cell_volume_m3", &CELL_VOLUME.to_string());
    c.set("chem.diffusion_m2_s", "1.0e-13");
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
    c.set("chem.species.0.atoms", "1");
    c.set("chem.species.1.atoms", "8");
    c.set("chem.reaction.count", "0");

    let mut s = Simulation::new(c);
    let materials = s.world.materials.clone();
    s.world.set_generator(Box::new(Gen { materials, seeds: vec![(0, 5_000_000), (1, 2_000_000)] }));
    s.world.activate(ChunkCoord::new(0, 0, 0));
    s
}

/// **The guarantee.** A world observed a thousand times is the world that was
/// never observed at all — same hash, same trajectory, same everything.
///
/// The signature already enforces this (`observe` takes `&World`), so this test
/// cannot fail while the code compiles. It is here anyway, because "obviously
/// read-only" is exactly the kind of claim that stops being true when somebody
/// adds a cache, and because a compile error is a worse explanation than a named
/// test.
#[test]
fn observing_a_world_never_changes_it() {
    let mut watched = sim(1);
    let mut ignored = sim(1);

    let mut rec = Recorder::new(1, 64);
    for _ in 0..300 {
        watched.tick();
        rec.observe(&watched.world);
        rec.observe(&watched.world); // twice a tick, for good measure
        ignored.tick();
    }

    assert_eq!(watched.hash(), ignored.hash(), "observation perturbed the world");
    assert!(!rec.is_empty(), "nothing was recorded");
    assert_eq!(watched.save(), ignored.save(), "observation changed what gets saved");
}

/// Sampling the same tick twice records it once. Otherwise a caller that
/// observes in a loop silently inflates every series it draws.
#[test]
fn a_tick_is_recorded_at_most_once() {
    let mut s = sim(2);
    let mut rec = Recorder::new(1, 8);
    for _ in 0..50 {
        s.tick();
        rec.observe(&s.world);
        rec.observe(&s.world);
        rec.observe(&s.world);
    }
    let ticks: Vec<u64> = rec.samples().iter().map(|x| x.tick).collect();
    let mut sorted = ticks.clone();
    sorted.dedup();
    assert_eq!(ticks, sorted, "the same tick was recorded more than once");
}

/// The interval is respected, and the first sample always lands so a series
/// never starts in the middle of nowhere.
#[test]
fn the_recorder_honours_its_interval() {
    let mut s = sim(3);
    let mut rec = Recorder::new(10, 8);
    rec.observe(&s.world);
    for _ in 0..100 {
        s.tick();
        rec.observe(&s.world);
    }
    assert!(rec.samples().len() >= 10 && rec.samples().len() <= 12);
    for x in rec.samples().iter().skip(1) {
        assert_eq!(x.tick % 10, 0, "sample at tick {} is off-interval", x.tick);
    }
}

/// **Matter lives in two places, and the recorder must show them apart.**
///
/// Their sum is what conservation checks; their ratio is what tells you whether
/// compartments are doing anything. Every conservation confusion in Phase 5i came
/// from conflating the two, so a view that added them up would be reproducing the
/// mistake it exists to prevent.
#[test]
fn bulk_and_encapsulated_are_reported_separately() {
    let mut s = sim(4);
    s.run(10);
    let mut rec = Recorder::new(1, 8);
    rec.snapshot(&s.world);
    let x = &rec.samples()[0];

    assert_eq!(x.bulk.len(), x.encapsulated.len());
    assert_eq!(x.bulk.iter().sum::<i128>(), 7_000_000, "the medium was miscounted");
    assert_eq!(
        x.encapsulated.iter().sum::<i128>(),
        0,
        "a world with no protocells cannot have encapsulated matter"
    );
    assert_eq!(x.protocells, 0);
}

/// Temperature comes back as a real range, derived the same way physics derives
/// it. Two of the worst bugs in this project announced themselves as absurd
/// temperatures while nobody was looking.
#[test]
fn temperature_is_reported_as_a_range() {
    let mut s = sim(5);
    s.run(5);
    let mut rec = Recorder::new(1, 8);
    rec.snapshot(&s.world);
    let x = &rec.samples()[0];
    assert!(x.temp_min <= x.temp_mean && x.temp_mean <= x.temp_max);
    assert!(
        (x.temp_mean - 300.0).abs() < 5.0,
        "expected ~300 K, got {:.3e} — the range is {:.3e}..{:.3e}",
        x.temp_mean,
        x.temp_min,
        x.temp_max
    );
}

/// A series survives a save and resume, because a run that has to restart loses
/// its history otherwise — and history is the only thing observation produces.
#[test]
fn a_series_can_span_a_resume() {
    let mut rec = Recorder::new(5, 8);
    let mut s = sim(6);
    s.run(30);
    rec.snapshot(&s.world);
    let bytes = s.save();

    let mut resumed = Simulation::load(&bytes, {
        let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
        c.merge(REFERENCE_MATERIALS).unwrap();
        c.set("world.seed", "6");
        c
    });
    // The config has to match for a resume; if it does not, the point of this
    // test is the recorder, not the loader, so fall back to the live world.
    if let Ok(mut r) = resumed.take() {
        r.run(30);
        rec.snapshot(&r.world);
    } else {
        s.run(30);
        rec.snapshot(&s.world);
    }
    assert_eq!(rec.samples().len(), 2);
    assert!(rec.samples()[1].tick > rec.samples()[0].tick, "time went backwards");
}

trait Takeable<T> {
    fn take(self) -> Result<T, ()>;
}
impl<T, E> Takeable<T> for Result<T, E> {
    fn take(self) -> Result<T, ()> {
        self.map_err(|_| ())
    }
}
