//! The open world — Phase 5h's integration.
//!
//! Every previous phase ran in a box. Energy could cross a boundary from Phase 2
//! onward, but matter never could: the reactor, the vent and the RAF survey all
//! took a fixed inventory of atoms and rearranged it. That is why the Phase 5c
//! survey was careful to claim only that its autocatalytic set *could* sustain
//! itself if fed. Nothing had ever fed anything.
//!
//! These tests are about the difference between the two things a system with a
//! gradient can settle into:
//!
//!   * **equilibrium** — the state a closed box reaches, where the gradient is
//!     gone, nothing flows, and nothing further happens, ever;
//!   * **a steady state** — the state an open system reaches, where the
//!     *contents* stop changing while matter continues to pour through, and the
//!     gradient is held up permanently by the flow.
//!
//! Only the second one can host anything alive. A dissipative structure is paid
//! for out of throughput, and a box has none. The analytic anchor is the twin of
//! the law Phase 2 validated conduction against: **steady one-dimensional
//! diffusion between a fixed source and a fixed sink gives a linear profile**
//! (Fick, 1855, in the same form Fourier gave heat in 1822), with a flux that is
//! the same across every face of the tube.

use au_core::hash::WorldHash;
use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, Lod};
use au_physics::{energy_for_temperature, Energy, GridSpec, Mass};
use au_sim::chem_columns::{install as install_chem, species_column};
use au_sim::physics_columns::{install as install_phys, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};
use au_sim::{Config, Simulation, DEFAULT_CONFIG};

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");

const NX: usize = 21;
const SOURCE_CELL: usize = 0;
const SINK_CELL: usize = NX - 1;
const HELD: i128 = 1_000_000;

/// A tube of molten host material, empty of species. Everything that ever
/// appears in it will have come through a port.
struct TubeGen {
    spec: GridSpec,
    pool: usize,
    temp_k: f64,
    materials: au_physics::MaterialRegistry,
}

impl ChunkGenerator for TubeGen {
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
        c
    }
}

/// Two elements, two monatomic species, no reactions — so the profile is pure
/// transport and the analytic result is not polluted by chemistry.
fn base_config(seed: u64) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("world.seed", &seed.to_string());
    c.set("clock.seconds_per_tick", "1.0");
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", &NX.to_string());
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", "1.0");
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("physics.gravity_m_s2", "0.0");
    c.set("chem.enabled", "true");
    c.set("chem.diffusion_m2_s", "0.05");
    c.set("chem.cell_volume_m3", "1.0");
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
    c
}

/// A tube with a source of species 0 at one end and an open port at the other.
fn open_config(seed: u64) -> Config {
    let mut c = base_config(seed);
    c.set("chem.reservoir.source", &format!("{}:0@{}", SOURCE_CELL, HELD));
    c.set("chem.reservoir.sink", &SINK_CELL.to_string());
    c
}

fn sim_with(config: Config) -> Simulation {
    let mut sim = Simulation::new(config);
    let spec = GridSpec::new(NX, 1, 1, 1.0);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(TubeGen { spec, pool: 2, temp_k: 2500.0, materials }));
    sim.world.chunk(ChunkCoord::new(0, 0, 0));
    sim.world.active.insert(ChunkCoord::new(0, 0, 0));
    sim
}

fn profile(sim: &Simulation, sp: usize) -> Vec<i128> {
    let chunk = sim.world.chunks.get(ChunkCoord::new(0, 0, 0)).unwrap();
    chunk.columns().get::<i128>(species_column(sp)).unwrap().as_slice().to_vec()
}

/// Atoms of each element currently in the world.
fn held(sim: &Simulation) -> Vec<i128> {
    let table_len = sim.world.atoms.len().max(2);
    let mut out = vec![0i128; table_len];
    for sp in 0..2usize {
        let m = sim.world.chem_registry.get(au_chem::SpeciesId(sp as u32)).unwrap();
        let total: i128 = profile(sim, sp).iter().sum();
        for (e, count) in m.formula(&elements()).into_iter().enumerate() {
            out[e] += count * total;
        }
    }
    out
}

fn elements() -> au_chem::PeriodicTable {
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
}

// ─── The books ───────────────────────────────────────────────────────────────

/// The identity the whole phase rests on, and the reason conservation does not
/// weaken when the world opens.
///
/// A sealed world asserts that totals never change. That was never the real
/// invariant — it was the special case of the real one where nothing crossed.
/// The real one is: **what is here now equals what was here at the start, plus
/// what came in, minus what went out** — exactly, as integers, per element.
///
/// If this can drift, something can eventually evolve to farm the drift. Free
/// matter is free food, and natural selection is an exhaustive search.
#[test]
fn the_books_balance_exactly_per_element() {
    let mut sim = sim_with(open_config(4));
    let start = held(&sim);
    for _ in 0..6_000 {
        sim.tick();
        let now = held(&sim);
        for e in 0..start.len() {
            assert_eq!(
                now[e] - start[e],
                sim.world.atoms.net(e),
                "element {}: the books disagree with the world",
                e
            );
        }
    }
    assert!(sim.world.atoms.inflow[0] > 0, "nothing was ever fed");
    assert!(sim.world.atoms.outflow[0] > 0, "nothing ever left");
}

/// A world with no ports books nothing, forever — so every earlier phase's
/// "totals are constant" reading still holds verbatim.
#[test]
fn a_sealed_world_never_books_a_single_atom() {
    let mut sim = sim_with(base_config(4));
    sim.run(500);
    assert!(!sim.world.atoms.any_flow());
    assert_eq!(sim.world.atoms.net(0), 0);
}

// ─── The boundary conditions themselves ──────────────────────────────────────

/// A source is a window onto something enormous: its cell stays at the declared
/// composition however fast the tube drains it.
#[test]
fn a_source_holds_its_cell_against_any_demand() {
    let mut sim = sim_with(open_config(11));
    for _ in 0..200 {
        sim.tick();
        assert_eq!(profile(&sim, 0)[SOURCE_CELL], HELD);
    }
}

/// A source tops up and never removes. Overfill its cell and it must leave the
/// excess alone — otherwise it is a sink in disguise and the books cannot tell
/// the two apart.
#[test]
fn a_source_never_removes_what_it_finds() {
    let mut sim = sim_with(open_config(12));
    {
        let chunk = sim.world.chunk(ChunkCoord::new(0, 0, 0));
        chunk.columns_mut().get_mut::<i128>(species_column(0)).unwrap().as_mut_slice()
            [SOURCE_CELL] = HELD * 10;
    }
    sim.tick();
    assert!(
        profile(&sim, 0)[SOURCE_CELL] > HELD,
        "the source drained a cell it was only supposed to fill"
    );
    assert_eq!(sim.world.atoms.outflow[0], 0, "a source booked an outflow");
}

/// **The constraint that matters most in this phase.**
///
/// A drain that removed some molecules and spared others would be a fitness
/// function written into a boundary condition — selection supplied as an input
/// rather than obtained as a result. So the sink is blind: whatever reaches it
/// is gone, regardless of what it is.
#[test]
fn a_sink_has_no_preferences() {
    // Transport off, so this is a statement about the port alone: hydrogen and
    // oxygen differ in mass by a factor of sixteen and diffuse at different
    // rates, and letting them move would confound "the drain does not choose"
    // with "the drain is reached at different speeds".
    let mut c = open_config(13);
    c.set("chem.diffusion_m2_s", "0.0");
    let mut open = sim_with(c);
    {
        let chunk = open.world.chunk(ChunkCoord::new(0, 0, 0));
        for sp in 0..2 {
            chunk.columns_mut().get_mut::<i128>(species_column(sp)).unwrap().as_mut_slice()
                [SINK_CELL] = 777 + sp as i128;
        }
    }
    open.tick();
    assert_eq!(profile(&open, 0)[SINK_CELL], 0);
    assert_eq!(profile(&open, 1)[SINK_CELL], 0);
    assert_eq!(open.world.atoms.outflow[0], 777, "hydrogen was treated differently");
    assert_eq!(open.world.atoms.outflow[1], 778, "oxygen was treated differently");
}

// ─── Steady state, and why it is not equilibrium ─────────────────────────────

/// The headline.
///
/// Closed: the tube flattens, the gradient dies, and nothing further happens.
/// Open: the contents stop changing *while matter keeps pouring through*, and
/// the profile that persists is the analytic one — **linear** between the held
/// end and the drained end, which is Fick's steady solution in one dimension.
///
/// The distinction is the whole reason this phase exists. Equilibrium is where a
/// box ends; a steady state is somewhere a system can *live*.
#[test]
fn an_open_tube_settles_into_a_linear_gradient_not_a_flat_equilibrium() {
    let mut sim = sim_with(open_config(21));
    sim.run(20_000);

    let p = profile(&sim, 0);
    // Interior cells only: the two port cells are boundary conditions, not
    // solutions of the equation.
    let n = SINK_CELL - SOURCE_CELL;
    let mut worst = 0.0f64;
    for i in SOURCE_CELL..=SINK_CELL {
        let expect = HELD as f64 * (1.0 - (i - SOURCE_CELL) as f64 / n as f64);
        let err = (p[i] as f64 - expect).abs() / HELD as f64;
        worst = worst.max(err);
    }
    assert!(worst < 0.02, "profile is not linear; worst deviation {:.4} of the held value", worst);

    // And it is a *flow*, not a rest state: the books keep moving while the
    // contents do not.
    let contents: i128 = p.iter().sum();
    let in_before = sim.world.atoms.inflow[0];
    sim.run(2_000);
    let after: i128 = profile(&sim, 0).iter().sum();
    let drift = (after - contents).abs() as f64 / contents as f64;
    assert!(drift < 0.01, "contents still changing by {:.4}; not a steady state", drift);
    assert!(
        sim.world.atoms.inflow[0] > in_before,
        "nothing flowed through — that is equilibrium, not a steady state"
    );
}

/// In a steady state the flux is the same across every face: whatever the source
/// supplies, the sink removes, with nothing accumulating anywhere in between.
#[test]
fn what_goes_in_comes_out() {
    let mut sim = sim_with(open_config(22));
    sim.run(20_000);
    let (i0, o0) = (sim.world.atoms.inflow[0], sim.world.atoms.outflow[0]);
    sim.run(5_000);
    let fed = sim.world.atoms.inflow[0] - i0;
    let lost = sim.world.atoms.outflow[0] - o0;
    assert!(fed > 0, "no throughput at all");
    let imbalance = (fed - lost).abs() as f64 / fed as f64;
    assert!(imbalance < 0.02, "fed {} but lost {} — that is accumulation", fed, lost);
}

// ─── Persistence and determinism ─────────────────────────────────────────────

/// An open world saved mid-flow and resumed must be the world that never
/// stopped — ledger included. Lose the ledger and a resumed world's books
/// silently disagree with its contents, which makes every later conservation
/// check a lie.
#[test]
fn an_open_world_survives_save_and_resume() {
    let mut straight = sim_with(open_config(33));
    straight.run(3_000);

    let mut resumed = sim_with(open_config(33));
    resumed.run(1_200);
    let bytes = resumed.save();
    let mut resumed = Simulation::load(&bytes, open_config(33)).unwrap();
    resumed.run(1_800);

    assert_eq!(straight.hash(), resumed.hash(), "a resumed open world is a different world");
    assert_eq!(straight.world.atoms, resumed.world.atoms, "the books did not survive the save");
    assert!(straight.world.atoms.any_flow());
}

/// Ports are part of the world's identity: two tubes differing only in whether
/// they are open must not agree on their history.
#[test]
fn opening_a_world_changes_what_it_is() {
    let mut sealed = sim_with(base_config(44));
    let mut open = sim_with(open_config(44));
    sealed.run(500);
    open.run(500);
    assert_ne!(sealed.hash(), open.hash());
}

// ─── Refusing malformed boundaries ───────────────────────────────────────────

/// A cell has one boundary condition or none. Two would fight, and which won
/// would depend on the order they were written in — a rule whose outcome depends
/// on its declaration order is not a rule.
#[test]
fn a_cell_cannot_be_both_a_source_and_a_sink() {
    let mut c = base_config(55);
    c.set("chem.reservoir.source", "3:0@100");
    c.set("chem.reservoir.sink", "3");
    let err = std::panic::catch_unwind(move || Simulation::new(c));
    assert!(err.is_err(), "a contradictory boundary was accepted");
}

/// A port naming a species the world has never heard of is a typo, and a typo in
/// a boundary condition produces a world that looks open, quietly runs to
/// equilibrium, and tells nobody why.
#[test]
fn a_port_cannot_hold_an_unknown_species() {
    let mut c = base_config(56);
    c.set("chem.reservoir.source", "0:99@100");
    let err = std::panic::catch_unwind(move || Simulation::new(c));
    assert!(err.is_err(), "a port named a species that does not exist");
}
