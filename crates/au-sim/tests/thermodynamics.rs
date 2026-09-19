//! Physics, inside a universe.
//!
//! `au-physics` is already proven against analytic solutions in its own crate.
//! What is tested here is everything that could go wrong in the *wiring* — and
//! the wiring is where the dangerous bugs are, because they are the ones that
//! look like nothing at all:
//!
//! * Does energy still balance to the microjoule once it is flowing through
//!   chunks, columns, snapshots and the scheduler?
//! * Does *looking* at the world still leave it unchanged, now that there is
//!   something in it that evolves?
//! * Does a saved world resume into the same future when the world has a
//!   thermal history?
//! * Does the clock still tell the truth, now that diffusion has an opinion
//!   about how fast it may run?

use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, Lod};
use au_physics::{energy_for_temperature, Energy, GridSpec, Mass, MaterialRegistry};
use au_sim::physics_columns::{install, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};
use au_sim::{Config, Simulation, DEFAULT_CONFIG};

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");

/// A slab of rock. An **initial condition**, not content — the same kind of thing
/// as the seed. Phase 4 will build planets; this builds a test rig.
struct SlabGen {
    spec: GridSpec,
    rho: f64,
    temp_k: f64,
    materials: MaterialRegistry,
}

impl ChunkGenerator for SlabGen {
    fn generate(&self, coord: ChunkCoord, _rng: &mut au_core::Rng) -> Chunk {
        let mut c = Chunk::new(coord, Lod::Full);
        let n = self.spec.cells();
        install(&mut c.columns, n);

        let id = self.materials.by_name("silicate").expect("fixture must define silicate");
        let mat = self.materials.get(id).unwrap();
        let kg = self.rho * self.spec.cell_volume();
        let e = energy_for_temperature(kg, self.temp_k, mat, 0.0);

        // Written directly, not through `columns_mut()`: generated state is
        // derivable from the seed and therefore not a *modification*. The chunk
        // stays clean and evictable — right up until physics touches it.
        c.columns.get_mut::<i128>(PHYS_MASS).unwrap().as_mut_slice().fill(Mass::from_kg(kg).0);
        c.columns
            .get_mut::<i128>(PHYS_ENERGY)
            .unwrap()
            .as_mut_slice()
            .fill(Energy::from_joules(e).0);
        c.columns.get_mut::<u16>(PHYS_MATERIAL).unwrap().as_mut_slice().fill(id.0);
        c
    }
}

fn slab_config(seed: u64, secs_per_tick: f64) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("world.seed", &seed.to_string());
    c.set("clock.seconds_per_tick", &secs_per_tick.to_string());
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", "2");
    c.set("physics.grid.ny", "16");
    c.set("physics.grid.nz", "2");
    c.set("physics.grid.cell_m", "10.0");
    c.set("physics.gravity_m_s2", "9.81");
    c.set("physics.bottom_flux_w_m2", "8.0");
    c.set("physics.top_radiates", "true");
    c
}

fn slab_sim(seed: u64, secs_per_tick: f64) -> Simulation {
    let cfg = slab_config(seed, secs_per_tick);
    let mut sim = Simulation::new(cfg);
    let spec = GridSpec::new(2, 16, 2, 10.0);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(SlabGen { spec, rho: 3000.0, temp_k: 300.0, materials }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));
    sim
}

fn total_energy(sim: &Simulation) -> i128 {
    sim.world
        .chunks
        .iter()
        .filter_map(|(_, c)| c.columns().get::<i128>(PHYS_ENERGY))
        .flat_map(|col| col.as_slice().iter())
        .sum()
}

// ═══ Conservation, end to end ════════════════════════════════════════════════

/// The whole point, restated at the level of a universe.
///
/// Energy flows in at the floor, out at the sky, and through twenty thousand
/// ticks of chunked, column-stored, scheduler-dispatched machinery. At the end,
/// the books balance **exactly** — not nearly, not to within a tolerance. Zero.
///
/// A single microjoule of discrepancy here is not an accuracy problem. It is a
/// free-energy exploit, and evolution is an exhaustive search that will find it.
#[test]
fn energy_balances_exactly_through_the_whole_engine() {
    let mut sim = slab_sim(7, 1.0e6);
    let e0 = total_energy(&sim);
    sim.run(20_000);
    let e1 = total_energy(&sim);

    let net = sim.world.ledger.net_energy().0;
    assert!(sim.world.ledger.energy_in.0 > 0, "no heat entered — the test is vacuous");
    assert!(sim.world.ledger.energy_out.0 > 0, "nothing radiated away");
    assert_eq!(
        e1 - e0,
        net,
        "ENERGY LEAK of {} µJ through the engine plumbing",
        (e1 - e0) - net
    );
}

/// Physics changes the world, so it must change the *hash* of the world. A
/// simulation whose fingerprint does not move is a simulation that is not
/// running.
#[test]
fn physics_actually_does_something() {
    let mut sim = slab_sim(7, 1.0e6);
    let h0 = sim.hash();
    let e0 = total_energy(&sim);
    sim.run(500);
    assert_ne!(sim.hash(), h0);
    assert_ne!(total_energy(&sim), e0, "no energy moved at all");
}

// ═══ The engine must still be invisible ══════════════════════════════════════

/// **The bug the active set exists to prevent.**
///
/// Phase 1 promised that looking at the world does not change it. Physics could
/// have destroyed that promise silently: if the physics system ran on "whatever
/// chunks are resident", then the set of things being *simulated* would be the
/// set of things somebody *looked at*, and the universe would evolve differently
/// depending on where the camera was pointed.
///
/// So physics runs on the **active set** — explicit, snapshotted, hashed — and
/// generating a chunk to look at does not enrol it. Here we generate fifty
/// chunks' worth of rock, in a world with physics running, and the history is
/// bit-identical.
#[test]
fn observing_the_world_does_not_simulate_it() {
    let mut quiet = slab_sim(21, 1.0e6);
    let mut watched = slab_sim(21, 1.0e6);

    // The camera sweeps across fifty chunks of real, physical rock. Each one is
    // generated. None is activated.
    for i in 1..51 {
        watched.world.chunk(ChunkCoord::new(i, 0, 0));
    }
    assert!(watched.world.chunks.resident() > quiet.world.chunks.resident());
    assert_eq!(watched.world.active.len(), 1, "looking must not enrol");

    quiet.run(3000);
    watched.run(3000);

    assert_eq!(
        quiet.hash(),
        watched.hash(),
        "the camera changed the universe"
    );
}

/// An active chunk is dirty forever — physics wrote to it, so it is no longer a
/// function of the seed and can never be regenerated. That is not a flaw; it is
/// the memory strategy telling the truth about what an evolving world costs.
///
/// It is also exactly why level of detail is not optional. Most of a planet must
/// *not* be finely evolving, or the planet will not fit in the machine.
#[test]
fn an_evolving_chunk_can_never_be_evicted_again() {
    let mut sim = slab_sim(3, 1.0e6);
    let coord = ChunkCoord::new(0, 0, 0);

    assert!(!sim.world.chunks.get(coord).unwrap().is_dirty(), "pristine before physics");
    sim.run(1);
    assert!(sim.world.chunks.get(coord).unwrap().is_dirty(), "physics must dirty what it touches");

    sim.world.evict_clean();
    assert!(sim.world.chunks.get(coord).is_some(), "a simulated chunk must survive eviction");
}

// ═══ Persistence ═════════════════════════════════════════════════════════════

/// A world with a thermal history must resume into the same future — including
/// its active set and its ledger. Lose the active set and a reloaded world
/// quietly stops simulating itself; lose the ledger and conservation can no
/// longer be checked across a save.
#[test]
fn a_thermal_history_survives_a_save() {
    let mut straight = slab_sim(99, 1.0e6);
    straight.run(6000);

    let mut split = slab_sim(99, 1.0e6);
    split.run(2500);
    let bytes = split.save();

    let mut resumed = Simulation::load(&bytes, slab_config(99, 1.0e6)).unwrap();
    // The generator is not persisted — it is a *rule*, not state, and the chunks
    // it made are in the snapshot. It is reinstalled only so that the resumed
    // world could grow, not to rebuild what it already has.
    let materials = resumed.world.materials.clone();
    resumed.world.set_generator(Box::new(SlabGen {
        spec: GridSpec::new(2, 16, 2, 10.0),
        rho: 3000.0,
        temp_k: 300.0,
        materials,
    }));
    resumed.run(3500);

    assert_eq!(resumed.world.active.len(), 1, "the active set must survive a save");
    assert_eq!(
        resumed.world.ledger, straight.world.ledger,
        "the books must survive a save"
    );
    assert_eq!(
        resumed.hash(),
        straight.hash(),
        "a reloaded universe took a different thermal path"
    );
}

// ═══ Time ════════════════════════════════════════════════════════════════════

/// **The finding that Phase 2 forced on Phase 1.**
///
/// The deep-time accelerator wants to double the seconds-per-tick forever.
/// Diffusion will not have it: an explicit heat solver diverges past
/// dt = dx²/(2·ndim·D), and no amount of wanting makes heat spread faster.
///
/// So the physics layer *caps the clock* and says so in the event log. The
/// accelerator keeps pushing; physics keeps refusing; the whole thing plateaus at
/// the truth.
///
/// The consequence is uncomfortable and worth stating plainly: **the deep-time
/// accelerator is largely inert on a world with real physics in it.** Reaching
/// evolutionary timescales will need coarse-grained physics at low LOD, not the
/// same physics run recklessly. Better to know in Phase 2 than in Phase 8.
#[test]
fn physics_refuses_to_let_the_clock_outrun_it() {
    // Ask for a hundred thousand years per tick. It cannot be done.
    let mut sim = slab_sim(5, 3.15e12);
    let asked = sim.world.clock.scale().as_secs_f64();
    sim.run(1);
    let allowed = sim.world.clock.scale().as_secs_f64();

    assert!(
        allowed < asked,
        "the clock ran at {:.3e} s/tick; diffusion cannot be integrated that fast",
        allowed
    );
    assert!(allowed > 0.0, "and it must not stop the clock entirely");

    let capped: Vec<_> = sim
        .world
        .events
        .by_layer(au_core::Layer::Physics)
        .filter(|e| e.kind == au_core::event::kind::PHYSICS_TIME_CAPPED)
        .collect();
    assert!(!capped.is_empty(), "the engine must SAY that it could not do what was asked");

    // And, having refused, it must still be stable. The point of the cap is that
    // the alternative is not a slightly-worse answer — it is temperatures
    // oscillating to infinity.
    sim.run(200);
    let e = total_energy(&sim);
    assert!(e > 0, "the grid went unstable and energy went negative");
}

/// An empty world imposes no speed limit — which is why every Phase 1 test still
/// passes untouched. Physics only slows the clock where there is physics to do.
#[test]
fn an_empty_world_still_accelerates_freely() {
    let mut sim = au_sim::boot(Some(1));
    sim.run(300_000);
    assert!(
        sim.world.clock.scale().secs > 1,
        "with nothing to simulate, there is nothing to slow down for"
    );
    assert!(sim.world.active.is_empty());
}

// ═══ Determinism, with physics running ═══════════════════════════════════════

#[test]
fn a_thermodynamic_universe_is_still_reproducible() {
    let trace = |seed: u64| {
        let mut s = slab_sim(seed, 1.0e6);
        (0..800).map(|_| { s.tick(); s.hash() }).collect::<Vec<_>>()
    };
    assert_eq!(trace(11), trace(11));
}

/// Structure, from nothing but a gradient.
///
/// Heat in at the floor, out at the sky. Nobody says "make the bottom hot". The
/// engine is not told there should be a temperature profile, or where it should
/// be, or what shape. It is told only that energy enters here and leaves there —
/// and a *layered world* appears, because that is what thermodynamics does with a
/// gradient.
///
/// This is small. It is also the first time in this project that something exists
/// which nobody put there, and every ambition in PROJECT_VISION is a bet that the
/// same trick keeps working as the layers pile up.
#[test]
fn a_gradient_produces_structure_that_nobody_designed() {
    let mut sim = slab_sim(1, 2.0e6);
    sim.run(30_000);

    let chunk = sim.world.chunks.get(ChunkCoord::new(0, 0, 0)).unwrap();
    let spec = GridSpec::new(2, 16, 2, 10.0);
    let e = chunk.columns().get::<i128>(PHYS_ENERGY).unwrap().as_slice();

    // Warmest at the floor, coldest at the sky, monotonic in between. Nobody
    // wrote that ordering down.
    for y in 1..spec.ny {
        let below = e[spec.idx(0, y - 1, 0)];
        let here = e[spec.idx(0, y, 0)];
        assert!(
            below > here,
            "the column is not stratified at y={} — no structure emerged",
            y
        );
    }
    let bottom = Energy(e[spec.idx(0, 0, 0)]).as_joules();
    let top = Energy(e[spec.idx(0, spec.ny - 1, 0)]).as_joules();
    assert!(bottom / top > 1.5, "the gradient is too weak to have organised anything");
}
