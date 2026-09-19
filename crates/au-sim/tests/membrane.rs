//! Membranes in the running world.
//!
//! au-chem proved the partition (free monomer buffered at the CMC) and the
//! geometry (a 1 µm³ bag needs ~10⁷ amphiphiles, which is what a real bacterium
//! has). au-physics proved the barrier (symmetric, exactly conservative). These
//! tests prove the tower: that a cell holding enough surfactant becomes hard to
//! reach, that this costs nothing, and that a membrane — being *derived* from
//! the populations rather than stored beside them — survives save and resume
//! without anyone having written it down.
//!
//! # The scale these tests run at
//!
//! A membrane means something at the scale where a plausible number of
//! amphiphiles can close a boundary. The geometry says that is a bacterium: a
//! 1 µm³ bag needs ~1.2×10⁷ molecules, and *E. coli* has about that many lipids.
//!
//! The engine could not represent that cell until the mass quantum was fixed.
//! Mass was stored as i128 *micrograms*, and a cubic micron of rock weighs
//! 2.5×10⁻¹⁵ kg — a millionth of one quantum. It rounded to zero, the host had
//! no mass, its phase could not be derived, the mobility mask called the cell
//! immobile, and nothing moved at all. The first version of these tests ran at
//! 1 µm and reported that "enabling membranes changed nothing"; what had really
//! happened was that not a single molecule had diffused anywhere.
//!
//! The quantum is now attograms, a bacterium-scale cell of rock is 2.5×10⁶ of
//! them, and these tests run where biology does.

use au_core::hash::WorldHash;
use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, Lod};
use au_physics::{energy_for_temperature, Energy, Mass};
use au_sim::chem_columns::{install as install_chem, species_column};
use au_sim::physics_columns::{install as install_phys, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};
use au_sim::{Config, Simulation, DEFAULT_CONFIG};

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");
const NX: usize = 9;
const CENTER: usize = 4;
const POOL: usize = 2;

/// A cube one micron on a side — the volume of a bacterium.
const CELL_M: f64 = 1.0e-6;
const CELL_VOLUME: f64 = 1.0e-18;
/// Roughly twice what it takes to close that volume (~1.2×10⁷ — which is also,
/// to an order of magnitude, the lipid count of a real bacterial membrane).
const AMPHIPHILE: i128 = 25_000_000;
const TRACER: i128 = 1_000_000;

/// Config strings for a straight-chain alcohol `CH₃(CH₂)ₙ₋₁OH`.
fn alcohol_spec(n: usize) -> (String, String) {
    let mut z: Vec<u32> = vec![6; n];
    let mut bonds: Vec<(usize, usize)> = (0..n.saturating_sub(1)).map(|i| (i, i + 1)).collect();
    let o = z.len();
    z.push(8);
    bonds.push((n - 1, o));
    let h = z.len();
    z.push(1);
    bonds.push((o, h));
    for i in 0..n {
        let deg = bonds.iter().filter(|(a, b)| *a == i || *b == i).count();
        for _ in deg..4 {
            let k = z.len();
            z.push(1);
            bonds.push((i, k));
        }
    }
    (
        z.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(","),
        bonds.iter().map(|(a, b)| format!("{}-{}:1", a, b)).collect::<Vec<_>>().join(","),
    )
}

struct TubeGen {
    temp_k: f64,
    seeds: Vec<(usize, i128)>,
    materials: au_physics::MaterialRegistry,
}

impl ChunkGenerator for TubeGen {
    fn generate(&self, coord: ChunkCoord, _rng: &mut au_core::Rng) -> Chunk {
        let mut c = Chunk::new(coord, Lod::Full);
        install_phys(&mut c.columns, NX);
        install_chem(&mut c.columns, NX, POOL);
        let id = self.materials.by_name("silicate").unwrap();
        let mat = self.materials.get(id).unwrap();
        let kg = 2.5e-15; // a cubic micron of rock: 2.5 × 10⁶ attograms
        let e = energy_for_temperature(kg, self.temp_k, mat, 0.0);
        c.columns.get_mut::<i128>(PHYS_MASS).unwrap().as_mut_slice().fill(Mass::from_kg(kg).0);
        c.columns
            .get_mut::<i128>(PHYS_ENERGY)
            .unwrap()
            .as_mut_slice()
            .fill(Energy::from_joules(e).0);
        c.columns.get_mut::<u16>(PHYS_MATERIAL).unwrap().as_mut_slice().fill(id.0);
        for &(sp, n) in &self.seeds {
            c.columns.get_mut::<i128>(species_column(sp)).unwrap().as_mut_slice()[CENTER] = n;
        }
        c
    }
}

/// A tube of molten silicate holding two declared species and no reactions: a
/// small tracer, and a five-carbon alcohol that is a genuine surfactant. Nothing
/// is created or destroyed — the only question is who can get out of the middle.
fn tube_config(membrane: bool, volume: f64) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("world.seed", "5041");
    c.set("clock.seconds_per_tick", "1e-5");
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", &NX.to_string());
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", &format!("{}", CELL_M));
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("physics.gravity_m_s2", "0.0");
    c.set("chem.enabled", "true");
    c.set("chem.cell_volume_m3", &format!("{}", volume));
    // D₀ of a 1 Da reference particle; a ~90 Da surfactant then moves at about
    // 10⁻⁹ m²/s, which is what a small molecule does in a liquid.
    c.set("chem.diffusion_m2_s", "1e-8");
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
    c.set("chem.species.count", "2");
    // 0: the tracer — hydroxyl. Two atoms, so no cut leaves a tail of two: it
    // can never assemble, whatever its concentration.
    c.set("chem.species.0.atoms", "1,8");
    c.set("chem.species.0.bonds", "0-1:1");
    // 1: pentanol — a genuine surfactant.
    let (atoms, bonds) = alcohol_spec(5);
    c.set("chem.species.1.atoms", &atoms);
    c.set("chem.species.1.bonds", &bonds);
    c.set("chem.reaction.count", "0");
    if membrane {
        c.set("chem.membrane", "true");
    }
    c
}

fn tube(membrane: bool, volume: f64, amphiphile: i128) -> Simulation {
    let mut sim = Simulation::new(tube_config(membrane, volume));
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(TubeGen {
        temp_k: 1500.0, // molten: the host is a fluid, so solutes may move
        seeds: vec![(0, TRACER), (1, amphiphile)],
        materials,
    }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));
    sim
}

fn profile(sim: &Simulation, sp: usize) -> Vec<i128> {
    let chunk = sim.world.chunks.iter().next().unwrap().1;
    chunk.columns().get::<i128>(species_column(sp)).unwrap().as_slice().to_vec()
}

/// **The headline.** Same tube, same tracer, same time — the only difference is
/// whether the surfactant in the middle cell is allowed to *be* a surface. When
/// it is, the tracer is still largely where it started.
///
/// Nothing decided that this cell should be a compartment. It holds a molecule
/// with a polar end and an apolar end, in enough quantity to exceed that
/// molecule's critical concentration and cover the boundary of a volume that
/// size, and everything else follows.
#[test]
fn a_cell_that_builds_a_membrane_keeps_its_contents() {
    let mut open = tube(false, CELL_VOLUME, AMPHIPHILE);
    let mut walled = tube(true, CELL_VOLUME, AMPHIPHILE);
    open.run(400);
    walled.run(400);

    let (o, w) = (profile(&open, 0), profile(&walled, 0));
    {
        let ch = open.world.chunks.iter().next().unwrap().1;
        let m = ch.columns().get::<i128>(au_sim::physics_columns::PHYS_MASS).unwrap().as_slice()[CENTER];
        let e = ch.columns().get::<i128>(au_sim::physics_columns::PHYS_ENERGY).unwrap().as_slice()[CENTER];
        }
    assert!(
        w[CENTER] > 3 * o[CENTER],
        "tracer retained: walled {} vs open {} — the membrane did nothing",
        w[CENTER],
        o[CENTER]
    );
    assert!(o[CENTER] < TRACER, "the open tube should have leaked");
    // The membrane holds its own material in, too — it is made of the thing it
    // is retaining, which is why a compartment that makes surfactant keeps it.
    let (ao, aw) = (profile(&open, 1), profile(&walled, 1));
    assert!(aw[CENTER] > ao[CENTER], "the surfactant should retain itself as well");
}

/// A surfactant too dilute to meet itself is just a solute. Spread the same
/// molecules through a far larger cell and they build nothing — and the world is
/// then bit-identical to one with membranes switched off entirely.
#[test]
fn a_dilute_amphiphile_builds_nothing() {
    let big = 1.0e-6; // the same molecules spread through a cubic centimetre
    let mut a = tube(false, big, AMPHIPHILE);
    let mut b = tube(true, big, AMPHIPHILE);
    a.run(200);
    b.run(200);
    assert_eq!(
        a.world.world_hash(),
        b.world.world_hash(),
        "below the critical concentration a membrane must not exist"
    );
}

/// Membranes are opt-in. An omitted flag and an explicit `false` give the same
/// world; turning it on changes it. Every earlier world is therefore untouched.
#[test]
fn membranes_are_off_by_default() {
    let mut a = tube(false, CELL_VOLUME, AMPHIPHILE);
    let mut b = Simulation::new({
        let mut c = tube_config(false, CELL_VOLUME);
        c.set("chem.membrane", "false");
        c
    });
    let materials = b.world.materials.clone();
    b.world.set_generator(Box::new(TubeGen {
        temp_k: 1500.0,
        seeds: vec![(0, TRACER), (1, AMPHIPHILE)],
        materials,
    }));
    b.world.activate(ChunkCoord::new(0, 0, 0));
    a.run(200);
    b.run(200);
    assert_eq!(a.world.world_hash(), b.world.world_hash(), "an explicit false differs from absent");

    let mut c = tube(true, CELL_VOLUME, AMPHIPHILE);
    c.run(200);
    assert_ne!(a.world.world_hash(), c.world.world_hash(), "enabling membranes changed nothing");
}

/// A barrier slows exchange; it never costs a molecule. Both species are
/// accounted for exactly across the whole tube.
#[test]
fn membranes_never_cost_a_molecule() {
    let mut sim = tube(true, CELL_VOLUME, AMPHIPHILE);
    sim.run(500);
    assert_eq!(profile(&sim, 0).iter().sum::<i128>(), TRACER, "the tracer leaked");
    assert_eq!(profile(&sim, 1).iter().sum::<i128>(), AMPHIPHILE, "the surfactant leaked");
    assert!(profile(&sim, 0).iter().all(|&v| v >= 0));
}

/// **Derived, never stored.** Because a membrane is a function of the
/// populations rather than a thing kept beside them, it needs no snapshot entry
/// and cannot fall out of step with the state it describes: a world saved
/// mid-diffusion and resumed rebuilds exactly the same barriers.
#[test]
fn a_membrane_survives_save_and_resume() {
    let mut sim = tube(true, CELL_VOLUME, AMPHIPHILE);
    sim.run(150);
    let bytes = sim.save();
    sim.run(150);
    let straight = sim.world.world_hash();

    let mut resumed =
        Simulation::load(&bytes, tube_config(true, CELL_VOLUME)).expect("resume");
    let materials = resumed.world.materials.clone();
    resumed.world.set_generator(Box::new(TubeGen {
        temp_k: 1500.0,
        seeds: vec![(0, TRACER), (1, AMPHIPHILE)],
        materials,
    }));
    resumed.run(150);
    assert_eq!(straight, resumed.world.world_hash(), "a membraned world diverged on resume");
}

/// Determinism is not negotiable, membranes or not.
#[test]
fn membraned_worlds_are_reproducible() {
    let run = || {
        let mut s = tube(true, CELL_VOLUME, AMPHIPHILE);
        s.run(200);
        s.world.world_hash()
    };
    assert_eq!(run(), run());
}

// ═══════════════════════════════════════════════════════════════════════════
//  Phase 5h — a membrane in a flow
// ═══════════════════════════════════════════════════════════════════════════

/// The same tube, now with an open port at one end drawing everything in its
/// cell to zero, forever.
fn draining_tube(membrane: bool, amphiphile: i128) -> Simulation {
    let mut c = tube_config(membrane, CELL_VOLUME);
    c.set("chem.reservoir.sink", &(NX - 1).to_string());
    let mut sim = Simulation::new(c);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(TubeGen {
        temp_k: 1500.0,
        seeds: vec![(0, TRACER), (1, amphiphile)],
        materials,
    }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));
    sim
}

/// **The point of the whole phase.**
///
/// A blind drain removes whatever reaches it, without preference — that is
/// enforced elsewhere and it is what keeps selection out of the boundary
/// condition. Put a compartment in front of one and something appears that
/// nobody wrote down: the cell that built a surface *keeps more of itself*.
///
/// Nothing in the engine grants membranes an advantage. Permeability multiplies
/// diffusive exchange (Phase 5g); a drain removes what reaches it (Phase 5h);
/// and a lower permeability therefore means less reaches the drain. The
/// advantage is what a barrier *is*, in the presence of a flow. It was already
/// implied by the two phases separately and became visible only when matter was
/// allowed to leave.
///
/// This is not yet selection — there is no heredity, no individuals, nothing
/// that could differentially reproduce. It is the precondition: an environment
/// in which two chemistries with identical composition have measurably different
/// prospects, decided by physics rather than by a parameter.
#[test]
fn a_membrane_is_an_advantage_only_once_the_world_is_open() {
    let mut walled = draining_tube(true, AMPHIPHILE);
    let mut open = draining_tube(false, AMPHIPHILE);
    walled.run(20_000);
    open.run(20_000);

    let w: i128 = profile(&walled, 0).iter().sum();
    let o: i128 = profile(&open, 0).iter().sum();
    assert!(
        w > 2 * o,
        "tracer surviving the drain: walled {} vs open {} — the membrane bought nothing",
        w,
        o
    );

    // And the difference is matter that *left*, not matter that was quietly
    // destroyed on one side and not the other. Checked on oxygen, because both
    // species carry exactly one: every molecule ever in either tube is one
    // oxygen atom, so the sum is a headcount.
    //
    // (The first draft of this checked hydrogen and compared the two worlds to
    // each other. Wrong twice: pentanol carries twelve hydrogens to the tracer's
    // one, and the two tubes are not supposed to hold the same amount at the end
    // — that is the entire result being tested.)
    for (name, sim) in [("walled", &walled), ("open", &open)] {
        let held: i128 =
            profile(sim, 0).iter().sum::<i128>() + profile(sim, 1).iter().sum::<i128>();
        assert_eq!(
            held + sim.world.atoms.outflow[2],
            TRACER + AMPHIPHILE,
            "{}: oxygen went missing",
            name
        );
    }
}

/// A dilute cell builds no surface, so it has nothing to hold anything back
/// with, and the drain treats it exactly as it treats a world with membranes
/// switched off. The advantage is a property of *having built* the barrier, not
/// of being permitted to.
#[test]
fn below_the_critical_concentration_a_membrane_buys_nothing_in_a_flow() {
    let mut dilute = draining_tube(true, AMPHIPHILE / 5_000);
    let mut none = draining_tube(false, AMPHIPHILE / 5_000);
    dilute.run(20_000);
    none.run(20_000);
    let d: i128 = profile(&dilute, 0).iter().sum();
    let n: i128 = profile(&none, 0).iter().sum();
    let diff = (d - n).abs() as f64 / n.max(1) as f64;
    assert!(diff < 0.05, "a dilute solution acted like a barrier: {} vs {}", d, n);
}

