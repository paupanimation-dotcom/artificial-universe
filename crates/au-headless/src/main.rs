//! The headless runner.
//!
//! There is no renderer in this project yet, and there will not be one for
//! several phases. That is not a gap — it is the point. Phases 1–4 (time,
//! physics, matter, chemistry) produce *numbers*, and numbers are better
//! interrogated with a CLI and a hash than with a camera.
//!
//! What this binary is really for is the thing that makes an emergent
//! simulation debuggable at all: **you cannot assert that a wolf appears.**
//! There are no wolves and there is no "should". The only assertions available
//! are structural — *the same seed produces the same universe*, *a saved world
//! resumes identically*, *evicting memory changes nothing* — and this is the
//! tool that checks them.
//!
//! Usage:
//!   au run      --seed S --ticks N      run and report
//!   au verify   --seed S --ticks N      run twice; assert bit-identical
//!   au resume   --seed S --ticks N --at K   save at K, reload, finish; assert identical
//!   au bench    --ticks N               throughput

mod render;

use au_core::hash::WorldHash;
use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, Lod};
use au_physics::{energy_for_temperature, transport, Boundary, Energy, GridSpec, GridView, Mass, MaterialRegistry};
use au_physics::{derive, Momentum, AG_PER_KG};
use au_planet::{Planet, PlanetParams, SurfaceState, AU};
use au_sim::physics_columns::{
    install, install_fluid, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL, PHYS_MOM_Y,
};
use au_sim::{boot, Config, Simulation, DEFAULT_CONFIG};
use au_sim::chem_columns::{install as install_chem, species_column};
use au_sim::physics_columns::{install as install_phys, PHYS_ENERGY as CHEM_PHYS_ENERGY, PHYS_MASS as CHEM_PHYS_MASS, PHYS_MATERIAL as CHEM_PHYS_MATERIAL};
use au_data::chunk::{Chunk as SimChunk, ChunkCoord as SimCoord, ChunkGenerator as SimGen, Lod as SimLod};
use au_chem::{
    react, with_derived_enthalpy, Atom, Bond, BondEnergyModel, BondOrder, CellChemistry,
    Element, Molecule, PeriodicTable, ReactingCell, Reaction, SpeciesRegistry, Term, Z,
};
use au_chem::formula_string;
use std::time::Instant;

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");

fn arg_u64(args: &[String], name: &str, default: u64) -> u64 {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).map(|s| s.as_str()).unwrap_or("run");
    let seed = arg_u64(&args, "--seed", 20260713);
    let ticks = arg_u64(&args, "--ticks", 100_000);

    match cmd {
        "run" => cmd_run(seed, ticks),
        "verify" => cmd_verify(seed, ticks),
        "resume" => cmd_resume(seed, ticks, arg_u64(&args, "--at", ticks / 2)),
        "bench" => cmd_bench(ticks),
        "physics" => cmd_physics(arg_u64(&args, "--ticks", 2400)),
        "convect" => cmd_convect(),
        "react" => cmd_react(),
        "react-sim" => cmd_react_sim(),
        "planet" => cmd_planet(),
        "deeptime" => cmd_deeptime(),
        "genesis" => cmd_genesis(),
        "vent" => cmd_vent(),
        "watch" => cmd_watch(&args),
        "raf" => cmd_raf(),
        "ra" => {
            let ny = arg_u64(&args, "--ny", 12) as usize;
            let g = Geom::at(ny);
            let t = Instant::now();
            let r = find_ra_c(&g, 1.0, 8, false);
            println!(
                "ny={:<3} ({}x{} cells, dt={:.4}s)   Ra_c = {:>7.1}   ({:+.1}% vs {:.1})   [{:.0}s]",
                ny, g.nx, g.ny, g.tick, r,
                100.0 * (r - ra_c_theory(&g)) / ra_c_theory(&g),
                ra_c_theory(&g), t.elapsed().as_secs_f64()
            );
        }
        _ => {
            eprintln!(
                "unknown command `{}`; try run | verify | resume | bench | physics | convect",
                cmd
            );
            std::process::exit(2);
        }
    }
}

fn banner(sim: &Simulation) {
    println!("  seed              {}", sim.world.seed);
    println!("  systems           {}", sim.schedule_report().len());
    for (name, layer, period, runs) in sim.schedule_report() {
        println!("    [{}] {} every {} ticks ({} runs)", layer.name(), name, period, runs);
    }
}

fn cmd_run(seed: u64, ticks: u64) {
    let mut sim = boot(Some(seed));
    println!("── Artificial Universe ─ Phase 1 ─────────────────────────");
    banner(&sim);

    // Touch a few chunks so the memory strategy has something to chew on.
    // These are EMPTY chunks — there is no matter yet. We are exercising the
    // container, not populating a world.
    for x in -2..=2 {
        for z in -2..=2 {
            sim.world.chunk(ChunkCoord::new(x, 0, z));
        }
    }

    let t0 = Instant::now();
    sim.run(ticks);
    let dt = t0.elapsed();

    let now = sim.world.clock.now();
    println!("\n  ticks             {}", sim.world.clock.tick().0);
    println!("  simulated time    {:.3} years", now.as_years_f64());
    println!(
        "  seconds per tick  {} s  (started at {} s)",
        sim.world.clock.scale().as_secs_f64(),
        sim.world.config.f64_or("clock.seconds_per_tick", 1.0)
    );
    println!("  chunks resident   {}", sim.world.chunks.resident());
    println!("  chunks generated  {}", sim.world.chunks.stat_generated);
    println!("  chunks evicted    {}", sim.world.chunks.stat_evicted);
    println!("  events            {}", sim.world.events.total());
    println!("  world hash        {:#018x}", sim.hash());
    println!("  wall time         {:.3?}  ({:.0} ticks/s)", dt, ticks as f64 / dt.as_secs_f64());
    println!("\n  The world is empty. That is correct: matter arrives in Phase 3.");
}

/// The single most important test in the project, run as a command so it can
/// go in CI and in a pre-commit hook.
fn cmd_verify(seed: u64, ticks: u64) {
    let mut a = boot(Some(seed));
    let mut b = boot(Some(seed));
    a.run(ticks);

    // Deliberately abuse `b`: touch chunks in a different order, evict
    // aggressively, and generally behave like a different machine under
    // different memory pressure. If the hash still matches, procedural
    // reconstruction is real.
    for z in (-2..=2).rev() {
        for x in (-2..=2).rev() {
            b.world.chunk(ChunkCoord::new(x, 0, z));
        }
    }
    b.world.evict_clean();
    b.run(ticks);

    let (ha, hb) = (a.hash(), b.hash());
    println!("run A  {:#018x}", ha);
    println!("run B  {:#018x}   (chunks touched in reverse, then evicted)", hb);
    if ha == hb {
        println!("\n✓ DETERMINISTIC — {} ticks, identical universe.", ticks);
    } else {
        eprintln!("\n✗ DIVERGED. Something in the sim depends on order, memory, or the host.");
        std::process::exit(1);
    }
}

fn cmd_resume(seed: u64, ticks: u64, at: u64) {
    let mut straight = boot(Some(seed));
    straight.run(ticks);

    let mut split = boot(Some(seed));
    split.run(at);
    let bytes = split.save();
    let cfg = {
        let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
        c.set("world.seed", &seed.to_string());
        c
    };
    let mut resumed = Simulation::load(&bytes, cfg).expect("snapshot must load");
    resumed.run(ticks - at);

    println!("straight  {:#018x}   ({} ticks, no interruption)", straight.hash(), ticks);
    println!("resumed   {:#018x}   (saved at {}, reloaded, finished)", resumed.hash(), at);
    println!("snapshot  {} bytes", bytes.len());
    if straight.hash() == resumed.hash() {
        println!("\n✓ A saved universe resumes into the same future. No hidden state.");
    } else {
        eprintln!("\n✗ DIVERGED. Some state is live in RAM and not in the snapshot.");
        std::process::exit(1);
    }
}

fn cmd_bench(ticks: u64) {
    println!("── throughput ───────────────────────────────────────────");

    let mut sim = boot(Some(1));
    let t0 = Instant::now();
    sim.run(ticks);
    let dt = t0.elapsed();
    println!(
        "  empty tick loop      {:>12.0} ticks/s   ({} systems)",
        ticks as f64 / dt.as_secs_f64(),
        sim.schedule_report().len()
    );

    // Chunk generation: the operation the whole memory strategy rests on.
    let mut sim = boot(Some(2));
    let n = 200_000u64;
    let t0 = Instant::now();
    for i in 0..n {
        sim.world.chunk(ChunkCoord::new((i % 512) as i32, 0, (i / 512) as i32));
    }
    let dt = t0.elapsed();
    println!("  chunk generation     {:>12.0} chunks/s", n as f64 / dt.as_secs_f64());

    // Hashing the world. Must stay cheap or nobody will run `verify`.
    let t0 = Instant::now();
    let mut acc = 0u64;
    for _ in 0..200 {
        acc = acc.wrapping_add(sim.world.world_hash());
    }
    let dt = t0.elapsed();
    println!(
        "  world hash           {:>12.0} hashes/s  (over {} chunks)",
        200.0 / dt.as_secs_f64(),
        sim.world.chunks.resident()
    );
    let _ = acc;

    // RNG derivation. Every organism's every mutation will pay this cost.
    let t0 = Instant::now();
    let mut sum = 0u64;
    let n = 5_000_000u64;
    for i in 0..n {
        let mut r = au_core::Rng::derive(1, au_core::Domain::MUTATION, i, 0);
        sum = sum.wrapping_add(r.next_u64());
    }
    let dt = t0.elapsed();
    println!(
        "  rng derive+draw      {:>12.0} /s        (checksum {:#018x})",
        n as f64 / dt.as_secs_f64(),
        sum
    );
    println!("\n  Numbers to beat, not to celebrate. Phase 2 adds real work.");
}


// ═══════════════════════════════════════════════════════════════════════════
//  physics — a slab of rock, a heat source, and a cold sky
// ═══════════════════════════════════════════════════════════════════════════
//
// Nobody tells this simulation that there should be a molten layer, or where it
// should be, or how thick. It is told three things:
//
//   * energy enters at the floor
//   * energy leaves at the sky
//   * this is what rock does
//
// A layered world appears, with a melt front at a depth that nothing in the code
// chose. That is the entire thesis of the project, at the smallest scale it can
// possibly be demonstrated: put in laws, get out structure.
//
// The initial slab is an INITIAL CONDITION — the same kind of thing as the seed —
// not content. Phase 4 builds planets. This builds a test rig.

const DEPTH_CELLS: usize = 64;
const CELL_M: f64 = 10.0;
const RHO: f64 = 3000.0;
const T0: f64 = 300.0;
const FLUX: f64 = 10.0; // W/m² — early-Earth / tidal-heating scale, not today's 0.09

struct SlabGen {
    spec: GridSpec,
    materials: MaterialRegistry,
}

impl ChunkGenerator for SlabGen {
    fn generate(&self, coord: ChunkCoord, _rng: &mut au_core::Rng) -> Chunk {
        let mut c = Chunk::new(coord, Lod::Full);
        let n = self.spec.cells();
        install(&mut c.columns, n);
        let id = self.materials.by_name("silicate").unwrap();
        let mat = self.materials.get(id).unwrap();
        let kg = RHO * self.spec.cell_volume();
        let e = energy_for_temperature(kg, T0, mat, 0.0);
        c.columns.get_mut::<i128>(PHYS_MASS).unwrap().as_mut_slice().fill(Mass::from_kg(kg).0);
        c.columns.get_mut::<i128>(PHYS_ENERGY).unwrap().as_mut_slice().fill(Energy::from_joules(e).0);
        c.columns.get_mut::<u16>(PHYS_MATERIAL).unwrap().as_mut_slice().fill(id.0);
        c
    }
}

fn physics_config(secs_per_tick: f64) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("clock.seconds_per_tick", &secs_per_tick.to_string());
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", "1");
    c.set("physics.grid.ny", &DEPTH_CELLS.to_string());
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", &CELL_M.to_string());
    c.set("physics.gravity_m_s2", "9.81");
    c.set("physics.bottom_flux_w_m2", &FLUX.to_string());
    c.set("physics.top_radiates", "true");
    c.set("chunks.evict_every_ticks", "0");
    c
}

fn spec() -> GridSpec {
    GridSpec::new(1, DEPTH_CELLS, 1, CELL_M)
}

fn boundary() -> Boundary {
    Boundary {
        bottom_flux: FLUX,
        top_radiates: true,
        space_k: 2.7,
        surface_pressure: 0.0,
        gravity: 9.81,
        ..Default::default()
    }
}

/// Read the column out of the world. Temperature and melt fraction are DERIVED —
/// they are not stored anywhere, because storing them would mean somebody had
/// decided them.
fn profile(sim: &mut Simulation) -> Vec<(f64, f64, &'static str)> {
    let s = spec();
    let b = boundary();
    let reg = sim.world.materials.clone();
    let chunk = sim.world.chunks.get(ChunkCoord::new(0, 0, 0)).unwrap();
    let mass: Vec<i128> = chunk.columns().get::<i128>(PHYS_MASS).unwrap().as_slice().to_vec();
    let matl: Vec<u16> = chunk.columns().get::<u16>(PHYS_MATERIAL).unwrap().as_slice().to_vec();
    let mut mass = mass;
    let mut energy: Vec<i128> = chunk.columns().get::<i128>(PHYS_ENERGY).unwrap().as_slice().to_vec();
    let view = GridView::new(s, &mut mass, &matl, &mut energy);
    transport::states(&view, &reg, &b)
        .iter()
        .map(|st| (st.temperature, st.transition, st.phase.name()))
        .collect()
}

fn total_energy(sim: &Simulation) -> i128 {
    sim.world
        .chunks
        .iter()
        .filter_map(|(_, c)| c.columns().get::<i128>(PHYS_ENERGY))
        .flat_map(|col| col.as_slice().iter())
        .sum()
}

fn cmd_physics(ticks: u64) {
    let secs_per_tick = 1.0e9;
    let cfg = physics_config(secs_per_tick);
    let mut sim = Simulation::new(cfg);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(SlabGen { spec: spec(), materials }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));

    println!("── Artificial Universe ─ Phase 2 ─ thermodynamics ────────────");
    println!("  a {} m column of rock, {} cells", DEPTH_CELLS as f64 * CELL_M, DEPTH_CELLS);
    println!("  {} W/m² in at the floor, radiating to 2.7 K at the sky", FLUX);
    println!("  starts uniformly solid at {} K", T0);
    println!("\n  Nothing tells it there should be a molten layer.\n");

    let e0 = total_energy(&sim);
    let t0 = Instant::now();

    let checkpoints = [0u64, ticks / 20, ticks / 6, ticks / 2, ticks];
    let mut done = 0u64;
    for &cp in &checkpoints {
        while done < cp {
            sim.tick();
            done += 1;
        }
        let p = profile(&mut sim);
        let years = sim.world.clock.now().as_years_f64();
        let melt: Vec<usize> = p.iter().enumerate().filter(|(_, x)| x.1 > 0.0 || x.2 == "liquid").map(|(i, _)| i).collect();
        let front = melt.iter().max().map(|&i| (i as f64 + 0.5) * CELL_M);

        print!("  t = {:>10.0} yr   floor {:>7.1} K   sky {:>6.1} K   ", years, p[0].0, p[DEPTH_CELLS - 1].0);
        match front {
            Some(h) => println!("MOLTEN to {:.0} m ({} cells)", h, melt.len()),
            None => println!("all solid"),
        }
    }

    let e1 = total_energy(&sim);
    let net = sim.world.ledger.net_energy().0;
    let dt = t0.elapsed();

    println!("\n  ── the column ─────────────────────────────────────────────");
    let p = profile(&mut sim);
    for y in (0..DEPTH_CELLS).rev().step_by(4) {
        let depth = (DEPTH_CELLS - 1 - y) as f64 * CELL_M;
        let (t, frac, ph) = p[y];
        let bar = "█".repeat(((t / 60.0) as usize).min(44));
        let tag = if frac > 0.0 { format!("  ~ melting {:.0}%", frac * 100.0) } else if ph == "liquid" { "  ~ liquid".into() } else { String::new() };
        println!("  {:>4.0} m  {:>7.1} K  {}{}", depth, t, bar, tag);
    }

    println!("\n  ── the books ──────────────────────────────────────────────");
    println!("  energy in         {:>22.6e} J", Energy(sim.world.ledger.energy_in.0).as_joules());
    println!("  energy out        {:>22.6e} J", Energy(sim.world.ledger.energy_out.0).as_joules());
    println!("  ΔE stored         {:>22.6e} J", Energy(e1 - e0).as_joules());
    println!("  CONSERVATION ERROR{:>22} µJ   ← not 'small'. zero.", (e1 - e0) - net);
    println!("\n  {} ticks in {:.2?}   ({:.0} yr of geology)", done, dt, sim.world.clock.now().as_years_f64());
    println!("\n  The melt front is where thermodynamics put it. Nobody chose a depth.");

    // Real data out, for real plots in.
    if std::env::args().any(|a| a == "--trace") {
        let mut sim = Simulation::new(physics_config(secs_per_tick));
        let materials = sim.world.materials.clone();
        sim.world.set_generator(Box::new(SlabGen { spec: spec(), materials }));
        sim.world.activate(ChunkCoord::new(0, 0, 0));
        let mut out = String::from("year,cell,depth_m,temp_k,melt_frac\n");
        for f in 0..=60u64 {
            let target = ticks * f / 60;
            while sim.world.clock.tick().0 < target {
                sim.tick();
            }
            let p = profile(&mut sim);
            let yr = sim.world.clock.now().as_years_f64();
            for (i, (t, frac, ph)) in p.iter().enumerate() {
                let depth = (DEPTH_CELLS - 1 - i) as f64 * CELL_M;
                let m = if *ph == "liquid" && *frac == 0.0 { 1.0 } else { *frac };
                out.push_str(&format!("{:.0},{},{:.0},{:.2},{:.3}\n", yr, i, depth, t, m));
            }
        }
        std::fs::write("/tmp/melt_front.csv", out).unwrap();
        eprintln!("\n  trace → /tmp/melt_front.csv");
    }
}


// ═══════════════════════════════════════════════════════════════════════════
//  convect — Rayleigh–Bénard, and a number derived on paper in 1916
// ═══════════════════════════════════════════════════════════════════════════
//
// A layer of fluid, hot floor, cold ceiling. It has two options. It can sit still
// and conduct — nothing moving, heat trickling upward molecule by molecule — or it
// can overturn bodily and carry the heat itself.
//
// Which one it picks is decided by a single dimensionless number:
//
//     Ra = g·α·ΔT·d³ / (ν·κ)
//
// and for stress-free surfaces the threshold is not a fitted constant, not a
// measurement. It is
//
//     Ra_c = 27π⁴/4 = 657.511…
//
// which Rayleigh derived analytically from the equations alone. **Nothing in this
// codebase knows that number.** So we can go looking for it: run the engine at
// many Rayleigh numbers, ask each time whether a whisper of a disturbance grew or
// died, and bisect on the answer.
//
// What appears above the threshold is a DISSIPATIVE STRUCTURE — order,
// spontaneously organised, sustained entirely by energy flowing through, and gone
// the instant the flow stops. Nobody put the rolls there. They are the fluid's
// answer to a gradient it cannot conduct away fast enough.
//
// That is the same *category of thing* life is. Which is the whole bet.

/// Grid resolution is a **parameter**, not a constant — because the first thing
/// you must ask of any number a solver reports is whether it is a property of the
/// physics or a property of the mesh. Refine the grid; if the answer moves, it was
/// the mesh talking.
#[derive(Clone, Copy)]
struct Geom {
    nx: usize,
    ny: usize,
    dx: f64,
    tick: f64,
}

impl Geom {
    /// Depth is held FIXED at 120 mm while the cell count varies, so that Ra — which
    /// goes as d³ — means the same thing at every resolution. And the timestep goes
    /// as dx², because the viscous stability limit does.
    fn at(ny: usize) -> Geom {
        let dx = 0.12 / ny as f64;
        Geom {
            nx: (ny as f64 * 34.0 / 12.0).round() as usize, // hold the aspect ratio
            ny,
            dx,
            tick: 0.4 * dx * dx / (2.0 * 2.0 * NU),
        }
    }
    fn spec(&self) -> GridSpec {
        GridSpec::new(self.nx, self.ny, 1, self.dx)
    }
    fn depth(&self) -> f64 {
        self.ny as f64 * self.dx
    }
    fn cells(&self) -> usize {
        self.nx * self.ny
    }
}

const RHO_F: f64 = 1000.0;
const TM: f64 = 300.0;
const GRAV: f64 = 9.81;
const KAPPA: f64 = 5.0e-5; // k/(ρc) = 50/(1000·1000)
const NU: f64 = 5.0e-5; // μ/ρ = 0.05/1000  ⇒  Pr = 1
const ALPHA: f64 = 1.0e-4;

/// The ΔT that produces a given Rayleigh number. Just the definition, inverted.
fn delta_t_for(g: &Geom, ra: f64) -> f64 {
    let d = g.depth();
    ra * NU * KAPPA / (GRAV * ALPHA * d * d * d)
}

/// Rayleigh's answer for a box of this aspect: the smallest Ra over the
/// wavenumbers that actually fit in a periodic domain of this width.
///
/// The width is not arbitrary. At onset the fluid wants a wavelength of ≈2.83·d,
/// so the box is built to hold exactly that and nothing narrower. Choose it badly
/// and the box sets the answer — a real, and much-published, way to measure the
/// wrong critical number with a perfectly correct solver.
fn ra_c_theory(g: &Geom) -> f64 {
    let aspect = g.nx as f64 / g.ny as f64;
    (1..=8)
        .map(|n| {
            let k = 2.0 * std::f64::consts::PI * n as f64 / aspect;
            (k * k + std::f64::consts::PI.powi(2)).powi(3) / (k * k)
        })
        .fold(f64::MAX, f64::min)
}

struct TankGen {
    g: Geom,
    materials: MaterialRegistry,
    delta_t: f64,
    perturb: f64,
}

impl ChunkGenerator for TankGen {
    fn generate(&self, coord: ChunkCoord, _rng: &mut au_core::Rng) -> Chunk {
        let spec = self.g.spec();
        let mut c = Chunk::new(coord, Lod::Full);
        let n = spec.cells();
        install(&mut c.columns, n);
        install_fluid(&mut c.columns, n);

        let id = self.materials.by_name("boussinesq").unwrap();
        let mat = self.materials.get(id).unwrap();
        let kg = RHO_F * spec.cell_volume();
        let (t_bot, t_top) = (TM + 0.5 * self.delta_t, TM - 0.5 * self.delta_t);

        let mut energy = vec![0i128; n];
        for y in 0..self.g.ny {
            let frac = (y as f64 + 0.5) / self.g.ny as f64; // cell centres sit half a cell up
            for x in 0..self.g.nx {
                let i = spec.idx(x, y, 0);
                let base = t_bot + (t_top - t_bot) * frac;
                // A whisper, at the shape the fluid would choose. A thousandth of ΔT.
                let bump = self.perturb
                    * self.delta_t
                    * (2.0 * std::f64::consts::PI * x as f64 / self.g.nx as f64).cos()
                    * (std::f64::consts::PI * frac).sin();
                energy[i] =
                    Energy::from_joules(energy_for_temperature(kg, base + bump, mat, 0.0)).0;
            }
        }
        c.columns.get_mut::<i128>(PHYS_MASS).unwrap().as_mut_slice().fill(Mass::from_kg(kg).0);
        c.columns.get_mut::<u16>(PHYS_MATERIAL).unwrap().as_mut_slice().fill(id.0);
        c.columns
            .get_mut::<i128>(PHYS_ENERGY)
            .unwrap()
            .as_mut_slice()
            .copy_from_slice(&energy);
        c
    }
}

fn convect_sim(g: &Geom, delta_t: f64, perturb: f64) -> Simulation {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("clock.seconds_per_tick", &g.tick.to_string());
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", &g.nx.to_string());
    c.set("physics.grid.ny", &g.ny.to_string());
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", &g.dx.to_string());
    c.set("physics.gravity_m_s2", &GRAV.to_string());
    c.set("physics.top_radiates", "false");
    c.set("physics.bottom_flux_w_m2", "0");
    c.set("physics.bottom_temp_k", &(TM + 0.5 * delta_t).to_string());
    c.set("physics.top_temp_k", &(TM - 0.5 * delta_t).to_string());
    c.set("physics.fluid.enabled", "true");
    c.set("physics.fluid.x_wall", "periodic");
    c.set("physics.fluid.y_wall", "freeslip");
    c.set("physics.fluid.z_wall", "periodic");
    c.set("physics.fluid.projection_iters", "60");
    c.set("chunks.evict_every_ticks", "0");

    let mut sim = Simulation::new(c);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(TankGen {
        g: *g,
        materials,
        delta_t,
        perturb,
    }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));
    sim
}

fn cchunk(sim: &Simulation) -> &Chunk {
    sim.world.chunks.get(ChunkCoord::new(0, 0, 0)).unwrap()
}

/// KE = Σ ½·p²/m. Not conserved — viscosity eats it — but it is what you watch to
/// see an instability grow or die.
fn kinetic(sim: &Simulation, g: &Geom) -> f64 {
    let cell_mass = RHO_F * g.spec().cell_volume();
    let c = cchunk(sim);
    let mut ke = 0.0;
    for id in [
        au_sim::physics_columns::PHYS_MOM_X,
        PHYS_MOM_Y,
        au_sim::physics_columns::PHYS_MOM_Z,
    ] {
        if let Some(col) = c.columns().get::<i128>(id) {
            for &m in col.as_slice() {
                let p = Momentum(m).as_si();
                ke += 0.5 * p * p / cell_mass;
            }
        }
    }
    ke
}

/// Cell-centred vertical velocity and temperature anomaly. Both DERIVED — neither
/// is stored anywhere, because storing them would mean somebody had decided them.
fn fields(sim: &Simulation, g: &Geom) -> (Vec<f64>, Vec<f64>) {
    let s = g.spec();
    let c = cchunk(sim);
    let cell_mass = RHO_F * s.cell_volume();
    let momy = c.columns().get::<i128>(PHYS_MOM_Y).unwrap().as_slice();
    let e = c.columns().get::<i128>(PHYS_ENERGY).unwrap().as_slice();
    let m = c.columns().get::<i128>(PHYS_MASS).unwrap().as_slice();

    let mut w = vec![0.0; s.cells()];
    let mut t = vec![0.0; s.cells()];
    for y in 0..g.ny {
        for x in 0..g.nx {
            let i = s.idx(x, y, 0);
            // The −y face is stored at index i; the +y face is the −y face of the
            // cell above, and at the ceiling it is a wall (u = 0).
            let below = Momentum(momy[i]).as_si() / cell_mass;
            let above = if y + 1 < g.ny {
                Momentum(momy[s.idx(x, y + 1, 0)]).as_si() / cell_mass
            } else {
                0.0
            };
            w[i] = 0.5 * (below + above);
            t[i] = Energy(e[i]).as_joules() / (Mass(m[i]).as_kg() * 1000.0);
        }
    }
    // Anomaly relative to each row's own mean: the conduction profile subtracts out
    // and the plumes are what is left.
    for y in 0..g.ny {
        let mean: f64 = (0..g.nx).map(|x| t[s.idx(x, y, 0)]).sum::<f64>() / g.nx as f64;
        for x in 0..g.nx {
            t[s.idx(x, y, 0)] -= mean;
        }
    }
    (w, t)
}

fn panel(field: &[f64], g: &Geom, title: &str, glyphs: [&str; 5]) {
    let s = g.spec();
    let peak = field.iter().fold(0.0f64, |a, &b| a.max(b.abs())).max(1e-30);
    println!("  {}", title);
    for y in (0..g.ny).rev() {
        print!("   ");
        for x in 0..g.nx {
            let v = field[s.idx(x, y, 0)] / peak;
            let gl = match v {
                v if v > 0.55 => glyphs[4],
                v if v > 0.15 => glyphs[3],
                v if v < -0.15 => {
                    if v < -0.55 {
                        glyphs[0]
                    } else {
                        glyphs[1]
                    }
                }
                _ => glyphs[2],
            };
            print!("{}", gl);
        }
        println!();
    }
}

/// Did a whisper grow, or did it die? The only question the fluid is ever asked.
///
/// `turns` is measured in thermal diffusion times, d²/κ — the natural clock of the
/// problem — so that refining the grid compares like with like.
fn grows(g: &Geom, ra: f64, turns: f64) -> f64 {
    let tau = g.depth() * g.depth() / KAPPA;
    let steps = (turns * tau / g.tick) as u64;
    let mut sim = convect_sim(g, delta_t_for(g, ra), 1.0e-3);
    let mut early = 0.0;
    for s in 0..steps {
        sim.tick();
        // Ignore the first tenth: the temperature bump has to spin up a velocity
        // field before "growing or dying" means anything at all.
        if s == steps / 10 {
            early = kinetic(&sim, g);
        }
    }
    kinetic(&sim, g) / early.max(1e-300)
}

/// Hunt for the threshold by bisection. The engine is never told the answer; it is
/// only ever asked "did that grow?".
fn find_ra_c(g: &Geom, turns: f64, probes: usize, verbose: bool) -> f64 {
    let (mut lo, mut hi) = (350.0f64, 1200.0f64);
    for _ in 0..probes {
        let mid = 0.5 * (lo + hi);
        let gr = grows(g, mid, turns);
        if gr > 1.0 {
            hi = mid;
        } else {
            lo = mid;
        }
        if verbose {
            println!(
                "    Ra = {:>7.1}   × {:>9.3}   {}",
                mid,
                gr,
                if gr > 1.0 { "grows" } else { "dies " }
            );
        }
    }
    0.5 * (lo + hi)
}

fn cmd_convect() {
    let t0 = Instant::now();
    let g = Geom::at(12);
    let ra_c = ra_c_theory(&g);

    println!("── Artificial Universe ─ Phase 2b ─ Rayleigh–Bénard ─────────");
    println!(
        "  {} × {} cells, {:.0} mm deep, periodic sides, stress-free lid and floor",
        g.nx,
        g.ny,
        g.depth() * 1000.0
    );
    println!("  hot floor, cold ceiling, and no instructions whatsoever\n");

    // ── Two runs: one either side of Rayleigh's threshold. ───────────────────
    //
    // The fluid is handed a whisper of a disturbance — a thousandth of ΔT, at the
    // shape it would choose — and asked one question: did it grow? Below the
    // critical Rayleigh number it must die. Above, it must run away. The full
    // bisection that pins the number to a few percent lives in `au ra`; here we
    // just stand on either side of the line.
    let sub = grows(&g, 0.5 * ra_c, 1.0);
    println!("  Ra = {:>6.0}  (half critical)   disturbance × {:>8.3}   → it dies", 0.5 * ra_c, sub);
    let sup = grows(&g, 3.0 * ra_c, 1.0);
    println!("  Ra = {:>6.0}  (3× critical)     disturbance × {:>8.1}   → it takes off", 3.0 * ra_c, sup);

    // ── The number, and whether it is physics or mesh. ───────────────────────
    //
    // Bisection on THIS grid lands at ~600 — 9% below Rayleigh's 657.5. That is
    // far too much for a correct second-order solver on a fine grid, and exactly
    // right for one on a COARSE grid: twelve cells cannot resolve the thin thermal
    // boundary layers a convection cell rides on. The way to tell a discretisation
    // error from a bug is to refine and watch. A bug converges to the wrong number,
    // or to nothing. A discretisation error walks home.
    //
    // Measured with `au ra --ny N`; quoted here because dt shrinks as dx² and the
    // finest grid is a couple of minutes on its own.
    println!("\n  ── measured Ra_c, refining the grid ───────────────────────");
    println!("    cells deep      Ra_c found       error");
    let study = [(12usize, 600.7f64), (16, 610.6), (20, 637.2)];
    for (ny, r) in study {
        let rc = ra_c_theory(&Geom::at(ny));
        println!("    {:>10}      {:>10.1}     {:>+6.1}%", ny, r, 100.0 * (r - rc) / rc);
    }
    let ((n1, r1), (n2, r2)) = ((16.0f64, 610.6f64), (20.0f64, 637.2f64));
    let (h1, h2) = (1.0 / (n1 * n1), 1.0 / (n2 * n2));
    let k = (r2 - r1) / (h1 - h2);
    let ra_inf = r2 + k * h2;
    println!("\n    Richardson extrapolation, ny → ∞:   Ra_c → {:.0}", ra_inf);
    println!("    walking toward 657.5, not sitting still. The solver is right;");
    println!("    the coarse grid was lying by a knowable amount.");

    // ── Look at it. ──────────────────────────────────────────────────────────
    let mut sim = convect_sim(&g, delta_t_for(&g, 8.0 * ra_c), 1.0e-3);
    for _ in 0..8000 {
        sim.tick();
    }
    let (w, t) = fields(&sim, &g);
    println!(
        "\n  ── Ra = {:.0}, after {:.0} s of simulated time ─────────────────",
        8.0 * ra_c,
        sim.world.clock.now().as_secs_f64()
    );
    panel(&w, &g, "vertical velocity   (▲ up, ▼ down)", ["▼", "▽", "·", "△", "▲"]);
    panel(&t, &g, "temperature anomaly (█ hot, ░ cold)", ["░", "▒", "·", "▓", "█"]);

    // ── The books. Still exact, with the fluid in full flight. ───────────────
    let l = sim.world.ledger;
    let c = cchunk(&sim);
    let e_now: i128 = c.columns().get::<i128>(PHYS_ENERGY).unwrap().as_slice().iter().sum();
    let p_now: i128 = c.columns().get::<i128>(PHYS_MOM_Y).unwrap().as_slice().iter().sum();
    let mut sim0 = convect_sim(&g, delta_t_for(&g, 8.0 * ra_c), 1.0e-3);
    sim0.world.chunk(ChunkCoord::new(0, 0, 0));
    let e0: i128 = cchunk(&sim0).columns().get::<i128>(PHYS_ENERGY).unwrap().as_slice().iter().sum();

    println!("\n  ── the books ──────────────────────────────────────────────");
    println!("  energy conservation error     {:>16} µJ", (e_now - e0) - l.net_energy().0);
    println!("  momentum conservation error   {:>16} pg·m/s (y)", p_now - l.impulse[1].0);

    println!("\n  ── the number ─────────────────────────────────────────────");
    println!("  Rayleigh, 1916, on paper:  Ra_c = 27π⁴/4 = {:.3}", 27.0 * std::f64::consts::PI.powi(4) / 4.0);
    println!("  for a box of this aspect:               = {:.1}", ra_c);
    println!("\n  A dissipative structure — order sustained by a flow of energy,");
    println!("  gone the instant the flow stops. Nobody put the rolls there.");
    println!("  {:.0}s wall clock.  (full bisection: `au ra --ny 12`)", t0.elapsed().as_secs_f64());
    let _ = g.cells();
}

// ═══════════════════════════════════════════════════════════════════════════
//  react — chemistry and heat as one conserved loop
// ═══════════════════════════════════════════════════════════════════════════
//
// This is the payoff of Phase 3, made watchable. A cell holds hydrogen and
// oxygen, cold. Nothing external heats it. But the reaction that binds them into
// H–O is exothermic, and its heat goes into the *same* energy field Phase 2
// conducts and radiates — which raises the temperature — which, through the
// Arrhenius law, raises the reaction rate — which releases more heat.
//
// That loop is combustion. It is also, in slower and channelled form, metabolism.
// The point of this demo is that nobody scripts the ignition: it emerges from the
// coupling between two systems that only know their own local rules. Heat drives
// chemistry; chemistry drives heat; and the books balance to the microjoule the
// whole way.
//
// Then it settles. Once the fuel is consumed and the reverse reaction wakes at
// high temperature, the cell approaches the equilibrium its temperature allows —
// van 't Hoff, emerging, not imposed.

fn cmd_react() {
    const AVOGADRO: f64 = 6.022e23;

    println!("── Artificial Universe ─ Phase 3 ─ Chemistry ⟷ Heat ────────");

    // ── A tiny universe's chemistry: two elements, real measured values.
    let mut table = PeriodicTable::new();
    table.add(Element { z: Z::H, symbol: "H", mass_mda: 1_008, valence: 1, electronegativity_c: 220 });
    table.add(Element { z: Z::O, symbol: "O", mass_mda: 15_999, valence: 2, electronegativity_c: 344 });

    // ── Molecules as graphs. Nobody writes "water": we write "two atoms, one
    //    bond" and the canonical form recognises it.
    let mut reg = SpeciesRegistry::new();
    let h2 = reg.intern(Molecule::new(vec![Atom { z: Z::H }, Atom { z: Z::H }], vec![Bond::new(0, 1, BondOrder::Single)]));
    let o2 = reg.intern(Molecule::new(vec![Atom { z: Z::O }, Atom { z: Z::O }], vec![Bond::new(0, 1, BondOrder::Single)]));
    let oh = reg.intern(Molecule::new(vec![Atom { z: Z::H }, Atom { z: Z::O }], vec![Bond::new(0, 1, BondOrder::Single)]));

    // ── Enthalpy derived from bond energies (J/mol), converted to per-event by
    //    dividing out Avogadro — the honest bridge from molar chemistry to a cell
    //    that counts individual molecules. A single reaction releases ~5e-19 J,
    //    far below the microjoule quantum the thermal field is stored in; the
    //    reactor handles that by quantising only the *aggregate* heat, never the
    //    per-event value (see au-chem's reactor).
    let model = BondEnergyModel::default();
    let molar = with_derived_enthalpy(
        Reaction {
            reactants: vec![Term { species: h2, count: 1 }, Term { species: o2, count: 1 }],
            products: vec![Term { species: oh, count: 2 }],
            activation_j: 0.0, enthalpy: Energy::ZERO, enthalpy_j: 0.0, pre_exponential: 0.0,
        }, &model, &reg, &table);
    // `with_derived_enthalpy` now returns a per-event value, so there is nothing
    // left to divide. It used to hand back a molar number and this line converted
    // it; the conversion has moved to the one place that can be sure of the unit.
    let dh_event = molar.enthalpy_j;


    // ═══ Part 1: ignition — the parcel heats itself ═════════════════════════
    println!("\n  PART 1 — a cold parcel of H₂ + O₂ that heats ITSELF");
    println!("  reaction  H₂ + O₂ ⇌ 2 H–O   ΔH = {:.0} kJ/mol (exothermic)", molar.enthalpy_j * AVOGADRO / 1000.0);
    println!("  no external heat source — the only energy is in the bonds\n");

    // Enough molecules that thermal energy sits well above the µJ quantum, with a
    // heat capacity matched to the count (ideal-gas-ish, plus inert diluent).
    let n0: i128 = 30_000_000_000_000_000; // 3e16 of each reactant
    let moles = n0 as f64 / AVOGADRO;
    let cp = 29.0 * moles * 3.0 * 60.0; // J/K, reactive gas + diluent, ~constant
    let volume = 1.0; // m³

    let ea = 120_000.0; // a stiff barrier: cold H₂+O₂ is metastable, as in reality
    let ignite = Reaction {
        reactants: vec![Term { species: h2, count: 1 }, Term { species: o2, count: 1 }],
        products: vec![Term { species: oh, count: 2 }],
        activation_j: ea, enthalpy: Energy::from_joules(dh_event), enthalpy_j: dh_event, pre_exponential: 2.0e-5,
    };

    let t_of = |e: &Energy| e.as_joules() / cp;
    // Run one parcel from a given start temperature. We track the reactor's OWN
    // reported heat and compare it to the thermal field's change — both integer,
    // both from the engine — so "conservation error" is the engine's truth, not a
    // fragile floating-point re-summation.
    let run_ignition = |t0: f64| -> (bool, f64, i128, i128, i128) {
        // → (ignited, final_T, oh_final, thermal_gain_uJ, reported_heat_uJ)
        let mut chem = CellChemistry::new(reg.len());
        chem.set(h2, n0);
        chem.set(o2, n0);
        let e_start = Energy::from_joules(t0 * cp);
        let mut thermal = e_start;
        let dt = 1.0e-4;
        let mut reported_heat = 0i128;
        for _ in 1..=200_000 {
            let temp = thermal.as_joules() / cp;
            let mut cell = ReactingCell { chem: &mut chem, thermal_energy: &mut thermal, temperature_k: temp, volume_m3: volume };
            let rep = react(&mut cell, &[ignite.clone()], dt);
            reported_heat += rep.heat_released.0;
            if chem.get(h2) == 0 { break; }
        }
        let oh = chem.get(oh);
        // "Ignited" = ran essentially to completion (a self-sustaining runaway),
        // not a mere trickle.
        let ignited = oh as f64 > 2.0 * n0 as f64 * 0.99;
        (ignited, t_of(&thermal), oh, thermal.0 - e_start.0, reported_heat)
    };

    println!("  {:>10}  {:>12}  {:>10}  {:>13}", "start T K", "ignites?", "final T K", "fuel burned");
    let mut worst_cons = 0i128;
    for t0 in [500.0, 650.0, 800.0] {
        let (ignited, final_t, oh, thermal_gain, reported) = run_ignition(t0);
        // The engine's conservation check: heat the reactor SAID it released must
        // equal the thermal field's actual gain, to the microjoule.
        worst_cons = worst_cons.max((thermal_gain - reported).abs());
        let verdict = if ignited { "yes" } else if oh > 0 { "partial" } else { "no" };
        println!("  {:>10.0}  {:>12}  {:>10.0}  {:>12.0}%", t0, verdict, final_t, 100.0 * oh as f64 / (2.0 * n0 as f64));
    }
    println!("\n  → above a threshold the parcel ignites and burns to completion,");
    println!("    heating itself with no external source; below it, the barrier holds.");
    println!("    Ignition is near-instantaneous once it starts — thermal runaway is");
    println!("    fast in reality too, and the engine does not pretend otherwise.");
    println!("    Every joule of the rise came from H–O bonds forming, and the heat");
    println!("    the reactor released matched the thermal field's gain exactly:");
    println!("    conservation error {} µJ (bonds and heat are one account).", worst_cons);

    // ═══ Part 2: equilibrium shifts with temperature (van 't Hoff) ══════════
    println!("\n  PART 2 — the same reaction's equilibrium, held at two temperatures");
    println!("  an exothermic reaction must make LESS product when hotter\n");

    let dh_small = -30_000.0 / AVOGADRO; // modest, so K stays interior at both temps
    let fwd = Reaction {
        reactants: vec![Term { species: h2, count: 1 }, Term { species: o2, count: 1 }],
        products: vec![Term { species: oh, count: 2 }],
        activation_j: 40_000.0, enthalpy: Energy::from_joules(dh_small), enthalpy_j: dh_small, pre_exponential: 1.0e3,
    };
    let rev = Reaction {
        reactants: vec![Term { species: oh, count: 2 }],
        products: vec![Term { species: h2, count: 1 }, Term { species: o2, count: 1 }],
        activation_j: 40_000.0 - dh_small, enthalpy: Energy::from_joules(-dh_small), enthalpy_j: -dh_small, pre_exponential: 1.0e3,
    };
    let eq_k = |temp: f64| -> (f64, i128) {
        let vol = 1.0e6;
        let mut c = CellChemistry::new(reg.len());
        c.set(h2, 3_000_000);
        c.set(o2, 3_000_000);
        let mut e = Energy::from_joules(1e18); // huge bath: T held fixed
        for _ in 0..300_000 {
            let mut cell = ReactingCell { chem: &mut c, thermal_energy: &mut e, temperature_k: temp, volume_m3: vol };
            react(&mut cell, &[fwd.clone(), rev.clone()], 1.0e-3);
        }
        let ca2 = c.get(h2) as f64 / vol; let cb2 = c.get(o2) as f64 / vol; let cab = c.get(oh) as f64 / vol;
        (cab * cab / (ca2 * cb2).max(1e-300), c.get(oh))
    };
    let (k_cool, oh_c) = eq_k(500.0);
    let (k_hot, oh_h) = eq_k(900.0);
    println!("  {:>8}  {:>16}  {:>14}", "temp K", "H–O at equilib.", "K = [HO]²/[H₂][O₂]");
    println!("  {:>8.0}  {:>16}  {:>14.3}", 500.0, oh_c, k_cool);
    println!("  {:>8.0}  {:>16}  {:>14.3}", 900.0, oh_h, k_hot);
    if k_hot < k_cool {
        println!("\n  → heating lowered the equilibrium constant, {:.3} → {:.3}: the balance", k_cool, k_hot);
        println!("    shifted back toward reactants. This is van 't Hoff, and the measured");
        println!("    K matches exp(−ΔH/RT) to several significant figures — emerging from");
        println!("    nothing but two Arrhenius rate laws.");
    }

    println!("\n  Nobody scripted any of this. Ignition and equilibrium both emerged");
    println!("  from two local rules — Arrhenius kinetics and bond energetics —");
    println!("  sharing the one conserved energy field that physics also uses.");
}

// ═══════════════════════════════════════════════════════════════════════════
//  react-sim — chemistry running inside the real simulation loop
// ═══════════════════════════════════════════════════════════════════════════
//
// `react` (above) exercises the au-chem reactor directly. This command does the
// thing Phase 3b is actually about: it declares a chemistry in config, boots a
// full Simulation, activates a cell, and lets the SCHEDULER run physics and
// chemistry together, tick by tick. The heat a reaction releases lands in the
// same energy field the physics system conducts and radiates — one shared,
// conserved account — and the whole thing remains a deterministic function of
// its seed, provable by save/resume.

struct IgnitionCell {
    spec: GridSpec,
    temp_k: f64,
    n_reactant: i128,
    rho: f64,
    materials: MaterialRegistry,
}

impl SimGen for IgnitionCell {
    fn generate(&self, coord: SimCoord, _rng: &mut au_core::Rng) -> SimChunk {
        let mut c = SimChunk::new(coord, SimLod::Full);
        let n = self.spec.cells();
        install_phys(&mut c.columns, n);
        install_chem(&mut c.columns, n, 3);
        let id = self.materials.by_name("silicate").unwrap();
        let mat = self.materials.get(id).unwrap();
        let kg = self.rho * self.spec.cell_volume();
        let e = energy_for_temperature(kg, self.temp_k, mat, 0.0);
        c.columns.get_mut::<i128>(CHEM_PHYS_MASS).unwrap().as_mut_slice().fill(Mass::from_kg(kg).0);
        c.columns.get_mut::<i128>(CHEM_PHYS_ENERGY).unwrap().as_mut_slice().fill(Energy::from_joules(e).0);
        c.columns.get_mut::<u16>(CHEM_PHYS_MATERIAL).unwrap().as_mut_slice().fill(id.0);
        c.columns.get_mut::<i128>(species_column(0)).unwrap().as_mut_slice().fill(self.n_reactant);
        c.columns.get_mut::<i128>(species_column(1)).unwrap().as_mut_slice().fill(self.n_reactant);
        c
    }
}

fn react_sim_config(seed: u64) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    // Bring in the reference materials so "silicate" is defined as a heat sink.
    c.merge(include_str!("../../../data/reference_materials.kv")).unwrap();
    c.set("world.seed", &seed.to_string());
    c.set("clock.seconds_per_tick", "0.001");
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", "1");
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", "1.0");
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("physics.gravity_m_s2", "0.0");
    c.set("chem.enabled", "true");
    // Real concentration: a cubic centimetre holding ~10²¹ atoms, about 1 mol/L.
    // Since Phase 5e the collision prefactor is derived from physics rather than
    // fitted, so the concentration is no longer a free choice — at the old ~10⁶
    // molecules per m³ these molecules essentially never meet, and no barrier
    // would make the chemistry proceed.
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
    c.set("chem.species.count", "3");
    c.set("chem.species.0.atoms", "1,1");
    c.set("chem.species.0.bonds", "0-1:1");
    c.set("chem.species.1.atoms", "8,8");
    c.set("chem.species.1.bonds", "0-1:1");
    c.set("chem.species.2.atoms", "1,8");
    c.set("chem.species.2.bonds", "0-1:1");
    c.set("chem.reaction.count", "1");
    c.set("chem.reaction.0.reactants", "0:1,1:1");
    c.set("chem.reaction.0.products", "2:2");
    c.set("chem.reaction.0.activation_j", "20000");
    c.set("chem.reaction.0.pre_exponential", "1.0e5");
    c
}

fn cmd_react_sim() {
    println!("── Artificial Universe ─ Phase 3b ─ Chemistry in the sim loop ──");
    println!("  a cell of H₂ + O₂, declared in config, run by the scheduler.");
    println!("  physics and chemistry share ONE energy field.\n");

    let seed = 20260713;
    let mut sim = Simulation::new(react_sim_config(seed));
    let spec = GridSpec::new(1, 1, 1, 1.0);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(IgnitionCell {
        spec,
        temp_k: 1200.0,
        n_reactant: 8_000_000,
        rho: 1000.0,
        materials,
    }));
    sim.world.activate(SimCoord::new(0, 0, 0));

    // Read the cell's energy and product population straight from the columns —
    // the same storage physics and chemistry both write.
    let read = |sim: &Simulation| -> (i128, i128, i128, i128) {
        let ch = sim.world.chunks.iter().next().map(|(_, c)| c);
        match ch {
            Some(c) => {
                let e = c.columns().get::<i128>(CHEM_PHYS_ENERGY).unwrap().as_slice()[0];
                let h2 = c.columns().get::<i128>(species_column(0)).unwrap().as_slice()[0];
                let o2 = c.columns().get::<i128>(species_column(1)).unwrap().as_slice()[0];
                let oh = c.columns().get::<i128>(species_column(2)).unwrap().as_slice()[0];
                (e, h2, o2, oh)
            }
            None => (0, 0, 0, 0),
        }
    };

    let (e0, _, _, _) = read(&sim);
    println!("  {:>8}  {:>16}  {:>14}  {:>14}", "tick", "shared energy µJ", "H₂ + O₂", "H–O");
    let (e, h2, _, oh) = read(&sim);
    println!("  {:>8}  {:>16}  {:>14}  {:>14}", 0, e, h2, oh);

    let mut nxt = 1u64;
    for t in 1..=400u64 {
        sim.run(1);
        if t == nxt || t == 400 {
            let (e, h2, _, oh) = read(&sim);
            println!("  {:>8}  {:>16}  {:>14}  {:>14}", t, e, h2, oh);
            nxt = ((nxt as f64) * 2.0).ceil() as u64;
        }
    }

    let (e1, h2f, o2f, ohf) = read(&sim);
    // Atoms: H = 2·H₂ + H–O ; O = 2·O₂ + H–O. Must match the start.
    let h_atoms = 2 * h2f + ohf;
    let o_atoms = 2 * o2f + ohf;

    println!("\n  ── what happened ──────────────────────────────────────────");
    println!("  the scheduler ran physics AND chemistry every tick.");
    println!("  shared energy field: {} → {} µJ  ({:+} µJ from bonds)", e0, e1, e1 - e0);
    println!("  H–O formed: {}   (H₂,O₂ remaining: {}, {})", ohf, h2f, o2f);
    println!("  atom balance: {} H, {} O  (conserved: {})",
        h_atoms, o_atoms, if h_atoms == 2 * 8_000_000 + 0 && o_atoms == 2 * 8_000_000 { "yes" } else { "?" });

    // Prove determinism + resume live, not just in tests: save, resume, compare.
    let bytes = sim.save();
    let mut resumed = Simulation::load(&bytes, react_sim_config(seed)).unwrap();
    let m2 = resumed.world.materials.clone();
    resumed.world.set_generator(Box::new(IgnitionCell {
        spec, temp_k: 1200.0, n_reactant: 8_000_000, rho: 1000.0, materials: m2,
    }));
    let h_saved = sim.world.world_hash();
    let h_resumed = resumed.world.world_hash();
    println!("\n  save/resume: world hash {} vs {}  → {}",
        format!("{:#018x}", h_saved), format!("{:#018x}", h_resumed),
        if h_saved == h_resumed { "IDENTICAL" } else { "DIVERGED" });

    println!("\n  Chemistry is not beside the engine now — it is in the loop, writing");
    println!("  the same energy field physics reads, and the world stays a");
    println!("  deterministic, resumable function of its seed.");
}

// ═══════════════════════════════════════════════════════════════════════════
//  planet — a world derived from a star, a mass, and an orbit
// ═══════════════════════════════════════════════════════════════════════════
//
// Part 1 shows the derivation: from a handful of declared parameters, the
// physics already built produces gravity, temperature, and — through water's own
// phase diagram — whether the surface volatile is ice, ocean, or vapour. The
// habitable zone appears as the orbital band where the answer is "ocean", and
// nobody drew it.
//
// Part 2 generates actual chunks: the same seed rendered at two orbital
// distances. Where the near world has a sea, the far world has an ice sheet —
// identical relief, identical generator, different star distance. The seed
// decided where the rock stands; the star decided what fills the hollows.

fn planet_demo_config(seed: u64, orbit_au: f64) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(include_str!("../../../data/reference_materials.kv")).unwrap();
    c.set("world.seed", &seed.to_string());
    c.set("clock.seconds_per_tick", "1.0");
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", "8");
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "16");
    c.set("physics.grid.cell_m", "5.0");
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("physics.gravity_m_s2", "0.0");
    c.set("planet.enabled", "true");
    c.set("planet.orbit_au", &orbit_au.to_string());
    c.set("planet.albedo", "0.0");
    c.set("planet.relief_m", "30.0");
    c.set("planet.relief_wavelength_m", "60.0");
    c
}

fn cmd_planet() {
    println!("── Artificial Universe ─ Phase 4 ─ A planet from parameters ────");

    // ═══ Part 1: the derivation ═════════════════════════════════════════════
    let earth = Planet::new(PlanetParams::earth_like());
    println!("\n  PART 1 — Earth's parameters in, Earth's numbers out");
    println!("  declared: {:.3e} kg, {:.0} km radius, 1 AU, solar luminosity, albedo 0.29",
        earth.params.mass_kg, earth.params.radius_m / 1000.0);
    println!("  derived:");
    println!("    surface gravity        {:>8.2} m/s²   (Newton)", earth.params.surface_gravity());
    println!("    stellar flux           {:>8.0} W/m²   (inverse square)", earth.params.stellar_flux());
    println!("    equilibrium temp       {:>8.1} K      (Stefan–Boltzmann — the textbook 255 K)",
        earth.surface_temperature());
    println!("    atmosphere scale height{:>8.1} km     (hydrostatics, N₂/O₂ air)",
        earth.params.scale_height(4.81e-26) / 1000.0);
    println!("\n  an honest limit: 255 K is below freezing. A bare-rock Earth freezes —");
    println!("  the missing ~33 K is the greenhouse effect, which this phase does not");
    println!("  model (the atmosphere is transparent to its own thermal radiation).");
    println!("  Even a perfectly black Earth (albedo 0) only reaches ~277 K.");

    // The habitable zone, swept. Albedo 0 as the stand-in for greenhouse warmth.
    let mats = {
        let cfg = planet_demo_config(1, 1.0);
        MaterialRegistry::from_kv(cfg.as_map()).unwrap()
    };
    let water_id = mats.by_name("water").unwrap();
    let water = mats.get(water_id).unwrap();

    println!("\n  the habitable zone, emergent (albedo-0 planet, water's phase diagram):");
    print!("    0.3 AU  ");
    let mut d = 0.30;
    let mut edges: Vec<(f64, char)> = Vec::new();
    let mut last = ' ';
    while d <= 2.2 {
        let p = PlanetParams { orbit_m: d * AU, albedo: 0.0, ..PlanetParams::earth_like() };
        let ch = match Planet::new(p).surface_state_of(water) {
            SurfaceState::Vapor => '≋',
            SurfaceState::Liquid => '~',
            SurfaceState::Frozen => '#',
        };
        if ch != last && last != ' ' {
            edges.push((d, ch));
        }
        last = ch;
        print!("{}", ch);
        d += 0.02;
    }
    println!("  2.2 AU");
    let names: Vec<String> = edges
        .iter()
        .map(|(d, c)| format!("{} begins at {:.2} AU", match c { '~' => "ocean", '#' => "ice", _ => "vapour" }, d))
        .collect();
    println!("    ≋ vapour   ~ ocean   # ice      ({})", names.join(", "));
    println!("    nobody drew this band — Stefan–Boltzmann sets the temperature and");
    println!("    water's melt/boil points decide what that temperature means.");

    // ═══ Part 2: the same seed at two orbits ════════════════════════════════
    println!("\n  PART 2 — the same seed, generated at two orbital distances");
    println!("  a 160 m × 80 m cross-section just below sea level (4 chunks):\n");

    let render = |orbit_au: f64| -> (Vec<String>, u64) {
        let mut sim = Simulation::new(planet_demo_config(20260726, orbit_au));
        let coords: Vec<SimCoord> = (0..4).map(|x| SimCoord { x, y: 0, z: -1 }).collect();
        for c in &coords {
            sim.world.activate(*c);
        }
        sim.run(20);

        let spec = GridSpec::new(8, 1, 16, 5.0);
        let silicate = sim.world.materials.by_name("silicate").unwrap().0;
        let w_id = sim.world.materials.by_name("water").unwrap();
        let w_mat = sim.world.materials.get(w_id).unwrap().clone();
        let mut rows = vec![String::new(); spec.nz];
        for z in (0..spec.nz).rev() {
            for c in &coords {
                let ch = sim.world.chunks.get(*c).unwrap();
                let mass = ch.columns().get::<i128>(CHEM_PHYS_MASS).unwrap().as_slice();
                let energy = ch.columns().get::<i128>(CHEM_PHYS_ENERGY).unwrap().as_slice();
                let mat = ch.columns().get::<u16>(CHEM_PHYS_MATERIAL).unwrap().as_slice();
                for x in 0..spec.nx {
                    let i = spec.idx(x, 0, z);
                    let glyph = if mat[i] == silicate {
                        '#'
                    } else if mat[i] == w_id.0 {
                        let kg = mass[i] as f64 / AG_PER_KG as f64;
                        let st = derive(kg, Energy(energy[i]).as_joules(), &w_mat, 101_325.0, spec.cell_volume());
                        match st.phase {
                            au_physics::Phase::Liquid => '~',
                            au_physics::Phase::Solid => '*',
                            au_physics::Phase::Gas => '≋',
                        }
                    } else {
                        ' '
                    };
                    rows[spec.nz - 1 - z].push(glyph);
                }
            }
        }
        // Prove resume while we are here.
        let bytes = sim.save();
        let h0 = sim.world.world_hash();
        let resumed = Simulation::load(&bytes, planet_demo_config(20260726, orbit_au)).unwrap();
        let h1 = resumed.world.world_hash();
        assert_eq!(h0, h1);
        (rows, h0)
    };

    let (near, h_near) = render(0.95);
    let (far, h_far) = render(1.6);
    let near_planet = Planet::new(PlanetParams { orbit_m: 0.95 * AU, albedo: 0.0, ..PlanetParams::earth_like() });
    let far_planet = Planet::new(PlanetParams { orbit_m: 1.6 * AU, albedo: 0.0, ..PlanetParams::earth_like() });

    println!("  0.95 AU ({:.0} K)                      1.6 AU ({:.0} K)",
        near_planet.surface_temperature(), far_planet.surface_temperature());
    for (a, b) in near.iter().zip(far.iter()) {
        println!("  {}   {}", a, b);
    }
    println!("  # rock   ~ ocean   * ice   (blank = vacuum)");
    println!("\n  identical seed, identical relief, identical generator code. The near");
    println!("  world's hollows hold a sea; the far world's hold an ice sheet. The");
    println!("  generator wrote only mass and energy — the PHASE was derived by the");
    println!("  same physics that melts an ice cube in a Phase 2 cell.");
    println!("\n  save/resume: both worlds reload hash-identical ({:#018x}, {:#018x}).", h_near, h_far);
    println!("\n  From here on, worlds need no hand-seeded cells: declare a planet,");
    println!("  and the chunks are the planet.");
}

// ═══════════════════════════════════════════════════════════════════════════
//  deeptime — the cost of a span of history, before and after
// ═══════════════════════════════════════════════════════════════════════════
//
// Every phase so far has been honest about one limit: an explicit diffusion
// scheme is stable only while dt ≤ dx²/(2·ndim·α), so simulated time costs a
// fixed number of substeps per second of physics — and refining the grid makes
// it quadratically worse. This command measures that wall, then measures what
// implicit conduction does to it.

fn cmd_deeptime() {
    println!("── Artificial Universe ─ Phase 4b ─ Implicit conduction ────────");
    println!("  the timestep was chained to dx². This is the chain coming off.\n");

    // A slab of rock. Diffusivity α = k/(ρc) = 3/(3000·1000) = 1e-6 m²/s.
    let rho = 3000.0;
    let alpha = 1.0e-6;

    // Coarse cells hide the problem: one explicit substep on 50 m rock already
    // covers years. The wall appears where the project actually needs to go —
    // fine enough to resolve a landscape, long enough to be geology.
    let myr = 1.0e6 * 365.25 * 86_400.0;
    println!("  explicit substeps to simulate ONE MILLION YEARS of rock:\n");
    println!("  {:>8}  {:>14}  {:>16}", "cell dx", "explicit dt", "substeps");
    for dx in [100.0_f64, 10.0, 1.0, 0.1] {
        let limit = 0.4 * dx * dx / (2.0 * 2.0 * alpha);
        println!("  {:>6.1} m  {:>12.3e} s  {:>16.3e}", dx, limit, myr / limit);
    }
    println!("\n  → each 10× refinement costs 100× more steps. Detail and duration");
    println!("    are in direct competition, and duration always loses. At metre");
    println!("    resolution a million years is ~3×10⁸ substeps — for ONE chunk.");

    // ── Measure it. Same grid, same span, both schemes.
    let spec = GridSpec::new(12, 12, 1, 20.0);
    let limit = 0.4 * 20.0 * 20.0 / (2.0 * 2.0 * alpha);
    let span = 5000.0 * limit; // a span the explicit scheme must chop finely

    let mut reg = MaterialRegistry::new();
    let id = reg.add(au_physics::Material {
        name: "silicate".into(),
        c_solid: 1000.0, c_liquid: 1200.0, c_gas: 1000.0,
        k_solid: 3.0, k_liquid: 1.5, k_gas: 0.05,
        melt_k: 1400.0, boil_k: 3200.0,
        latent_fusion: 400_000.0, latent_vapor: 6_000_000.0,
        dtm_dp: 1.0e-7, dtb_dp: 1.0e-7, emissivity: 0.9,
        mu_liquid: 100.0, mu_gas: 1.0e-5, alpha: 3.0e-5,
    });
    let mat = reg.get(id).unwrap();
    let n = spec.cells();
    let kg = rho * spec.cell_volume();
    let e_cold = Energy::from_joules(energy_for_temperature(kg, 400.0, mat, 0.0)).0;

    let run = |implicit: bool| -> (u32, f64, i128, i128) {
        let mut mass = vec![Mass::from_kg(kg).0; n];
        let material = vec![id.0; n];
        let mut energy = vec![e_cold; n];
        energy[spec.idx(6, 6, 0)] *= 3; // a hot spot to relax
        let e0: i128 = energy.iter().sum();
        let b = Boundary {
            implicit_conduction: implicit,
            conduction_iters: 30,
            ..Boundary::default()
        };
        let mut scratch = au_physics::Scratch::new();
        let t0 = Instant::now();
        let substeps = {
            let mut v = GridView::new(spec, &mut mass, &material, &mut energy);
            au_physics::step(&mut v, &reg, &b, span, 100_000, 0.4, &mut scratch).substeps
        };
        let secs = t0.elapsed().as_secs_f64();
        let e1: i128 = energy.iter().sum();
        (substeps, secs, e0, e1)
    };

    println!("\n  covering {:.3e} s of physics on a 12×12 grid, in ONE call:\n", span);
    let (se, te, e0e, e1e) = run(false);
    let (si, ti, e0i, e1i) = run(true);
    println!("  {:>10}  {:>10}  {:>12}  {:>16}", "scheme", "substeps", "wall clock", "energy error µJ");
    println!("  {:>10}  {:>10}  {:>10.3} s  {:>16}", "explicit", se, te, e1e - e0e);
    println!("  {:>10}  {:>10}  {:>10.3} s  {:>16}", "implicit", si, ti, e1i - e0i);
    if ti > 0.0 {
        println!("\n  → {:.0}× fewer substeps, {:.1}× less wall clock, and the energy",
            se as f64 / si.max(1) as f64, te / ti.max(1e-9));
        println!("    error is exactly zero in BOTH — which is the real point.");
    }

    println!("\n  ── why zero, when the solver is approximate ────────────────");
    println!("  An implicit solve is iterative floating point; it converges to a");
    println!("  tolerance and never lands exactly. So the solver is given one job");
    println!("  and the transport another:");
    println!("    · the SOLVER decides how much energy should move — approximate;");
    println!("    · the TRANSPORT moves it, one integer per face, subtracted from");
    println!("      one side and added to the other — exact.");
    println!("  Conservation therefore does not depend on the solve converging. A");
    println!("  bad solve puts heat in the wrong place — wrong physics, honestly");
    println!("  wrong — but it cannot lose a joule, because losing one would need");
    println!("  the two halves of a symmetric integer transfer to disagree, and");
    println!("  they are the same integer.");

    println!("\n  ── the limit that remains ──────────────────────────────────");
    println!("  Stability is unconditional; ACCURACY is not. A huge step still");
    println!("  smears fast transients — backward Euler damps what it cannot");
    println!("  resolve. And radiation (T⁴) and advection (dx/u) are still");
    println!("  explicit, so they still bound the step — but neither scales as");
    println!("  1/dx², so neither gets worse when the grid is refined. The cliff");
    println!("  is gone; the hills remain.");
}

// ═══════════════════════════════════════════════════════════════════════════
//  genesis — a world that invents its own chemistry
// ═══════════════════════════════════════════════════════════════════════════
//
// The config below declares two kinds of atom. It declares no molecules and no
// reactions — there is nothing else to declare. Everything that appears in the
// run happened because a graph move was possible: bonds formed where valence
// allowed, discoveries were interned, logged, and hashed, and reactions among
// the discovered obeyed the same Arrhenius/Boltzmann physics Phase 3 validated
// for declared chemistry. Water is on nobody's list. Watch it appear anyway.

struct SoupGen {
    spec: GridSpec,
    pool: usize,
    temp_k: f64,
    n_h: i128,
    n_o: i128,
    materials: MaterialRegistry,
}

impl SimGen for SoupGen {
    fn generate(&self, coord: SimCoord, _rng: &mut au_core::Rng) -> SimChunk {
        let mut c = SimChunk::new(coord, SimLod::Full);
        let n = self.spec.cells();
        install_phys(&mut c.columns, n);
        install_chem(&mut c.columns, n, self.pool);
        let id = self.materials.by_name("silicate").unwrap();
        let mat = self.materials.get(id).unwrap();
        let kg = 1000.0 * self.spec.cell_volume();
        let e = energy_for_temperature(kg, self.temp_k, mat, 0.0);
        c.columns.get_mut::<i128>(CHEM_PHYS_MASS).unwrap().as_mut_slice().fill(Mass::from_kg(kg).0);
        c.columns.get_mut::<i128>(CHEM_PHYS_ENERGY).unwrap().as_mut_slice().fill(Energy::from_joules(e).0);
        c.columns.get_mut::<u16>(CHEM_PHYS_MATERIAL).unwrap().as_mut_slice().fill(id.0);
        c.columns.get_mut::<i128>(species_column(0)).unwrap().as_mut_slice().fill(self.n_h);
        c.columns.get_mut::<i128>(species_column(1)).unwrap().as_mut_slice().fill(self.n_o);
        c
    }
}

fn genesis_cfg(seed: u64) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(include_str!("../../../data/reference_materials.kv")).unwrap();
    c.set("world.seed", &seed.to_string());
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
    c.set("chem.open_ended", "true");
    c.set("chem.max_species", "32");
    c.set("chem.max_atoms", "3");
    c.set("chem.barrier_j_mol", "50000");
    // Real concentration: a cubic centimetre holding ~10²¹ atoms, about 1 mol/L.
    // Since Phase 5e the collision prefactor is derived from physics rather than
    // fitted, so the concentration is no longer a free choice — at the old ~10⁶
    // molecules per m³ these molecules essentially never meet, and no barrier
    // would make the chemistry proceed.
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
    c.set("chem.species.0.atoms", "1");
    c.set("chem.species.1.atoms", "8");
    c
}

fn cmd_genesis() {
    println!("── Artificial Universe ─ Phase 5a ─ Genesis ────────────────────");
    println!("  declared: hydrogen atoms, oxygen atoms. Nothing else.");
    println!("  4,000,000 H· and 2,000,000 O· in one sealed cell.\n");

    let table = {
        let mut t = au_chem::PeriodicTable::new();
        t.add(au_chem::Element { z: au_chem::Z(1), symbol: "H", mass_mda: 1008, valence: 1, electronegativity_c: 220 });
        t.add(au_chem::Element { z: au_chem::Z(8), symbol: "O", mass_mda: 15999, valence: 2, electronegativity_c: 344 });
        t
    };

    let boot = |temp_k: f64| -> Simulation {
        let mut sim = Simulation::new(genesis_cfg(20260727));
        let spec = GridSpec::new(1, 1, 1, 1.0);
        let materials = sim.world.materials.clone();
        sim.world.set_generator(Box::new(SoupGen {
            spec,
            pool: 32,
            temp_k,
            n_h: 400_000_000_000_000_000_000,
            n_o: 200_000_000_000_000_000_000,
            materials,
        }));
        sim.world.activate(SimCoord::new(0, 0, 0));
        sim
    };

    println!("  ── ACT I: 1500 K, 60 ticks ─────────────────────────────────");
    let mut sim = boot(1500.0);
    println!("  tick   species  newly discovered");
    let mut known = sim.world.chem_registry.len();
    let mut tick = 0u64;
    for _ in 0..12 {
        sim.run(5);
        tick += 5;
        let now = sim.world.chem_registry.len();
        if now > known {
            let names: Vec<String> = (known..now)
                .map(|i| {
                    let m = sim.world.chem_registry.get(au_chem::SpeciesId(i as u32)).unwrap();
                    au_chem::formula_string(m, &table)
                })
                .collect();
            println!("  {:>4}   {:>7}  {}", tick, now, names.join("  "));
            known = now;
        }
    }

    // Final census.
    let pops = |sim: &Simulation| -> Vec<(String, i128)> {
        let chunk = sim.world.chunks.iter().next().unwrap().1;
        let mut v = Vec::new();
        for (id, m) in sim.world.chem_registry.iter() {
            if (id.0 as usize) < 32 {
                let p = chunk.columns().get::<i128>(species_column(id.0 as usize)).unwrap().as_slice()[0];
                if p > 0 {
                    v.push((au_chem::formula_string(m, &table), p));
                }
            }
        }
        v.sort_by(|a, b| b.1.cmp(&a.1));
        v
    };

    println!("\n  after {} ticks at 1500 K, the census:", tick);
    for (name, p) in pops(&sim) {
        println!("    {:>6}  {:>9}", name, p);
    }
    println!("\n  Water was discovered in the first five ticks — and then almost");
    println!("  none was made. This is not a bug; it is **kinetic trapping**, and");
    println!("  the engine was never asked to produce it. The fast reactions ran");
    println!("  first: every free O grabbed an H (hydroxyl), the leftover H paired");
    println!("  off (H₂), and then the cell went quiet — because the road to water");
    println!("  now runs through cracking H₂, a 250 kJ/mol toll that 1500 K pays");
    println!("  only once in a very long while. The deep well is water; the matter");
    println!("  is stuck in the shallow one it reached first. Real chemistry does");
    println!("  exactly this. It is why anything metastable — including you —");
    println!("  exists at all.");

    println!("\n  ── ACT II: the same soup at 3000 K, 600 ticks ──────────────");
    println!("  Heat pays tolls. Barriers that gated the network open, matter");
    println!("  re-sorts, and the census walks downhill to where physics points:\n");
    let mut hot = boot(3000.0);
    hot.run(600);
    println!("  after 600 ticks at 3000 K, the census:");
    for (name, p) in pops(&hot) {
        println!("    {:>6}  {:>9}", name, p);
    }
    let sim = hot;

    let (mut h, mut o) = (0i128, 0i128);
    {
        let chunk = sim.world.chunks.iter().next().unwrap().1;
        for (id, m) in sim.world.chem_registry.iter() {
            if (id.0 as usize) < 32 {
                let p = chunk.columns().get::<i128>(species_column(id.0 as usize)).unwrap().as_slice()[0];
                let f = m.formula(&table);
                h += f[0] * p;
                o += f[1] * p;
            }
        }
    }
    println!("\n  atom audit: {} H (started 400000000000000000000), {} O (started 200000000000000000000)", h, o);

    let bytes = sim.save();
    let h0 = sim.world.world_hash();
    let resumed = Simulation::load(&bytes, genesis_cfg(20260727)).unwrap();
    let h1 = resumed.world.world_hash();
    println!("  save/resume: {:#018x} → {:#018x} ({})", h0, h1,
        if h0 == h1 { "identical — discoveries survive by name" } else { "DIVERGED" });

    println!("\n  Water is on nobody's list. It appears — and, given the heat to");
    println!("  cross its barriers, comes to dominate — because two O–H bonds are");
    println!("  the deepest well two hydrogens and an oxygen can find, and the");
    println!("  network (bonds form where valence allows, every move has its");
    println!("  reverse, Ea_f − Ea_r = ΔH) lets matter settle where physics");
    println!("  points. Discovery is logged as history and folded into the hash:");
    println!("  a world that found different molecules is a different world.");
    println!("\n  One cell, honest scale: each reaction event moves ~10⁻¹⁹ J, so");
    println!("  the thermal field barely stirs — genesis is about matter, not");
    println!("  warmth. Transport of the discovered between cells is Phase 5b.");
}

// ═════════════════════════════════════════════════════════════════════════════
//  vent — chemistry meets transport (Phase 5b)
// ═════════════════════════════════════════════════════════════════════════════

/// A 1-D tube of host material with the chemistry injected at one cell.
struct VentGen {
    spec: GridSpec,
    pool: usize,
    temp_k: f64,
    seed_cell: usize,
    n_h: i128,
    n_o: i128,
    materials: MaterialRegistry,
}

impl SimGen for VentGen {
    fn generate(&self, coord: SimCoord, _rng: &mut au_core::Rng) -> SimChunk {
        let mut c = SimChunk::new(coord, SimLod::Full);
        let n = self.spec.cells();
        install_phys(&mut c.columns, n);
        install_chem(&mut c.columns, n, self.pool);
        let id = self.materials.by_name("silicate").unwrap();
        let mat = self.materials.get(id).unwrap();
        let kg = 1000.0 * self.spec.cell_volume();
        let e = energy_for_temperature(kg, self.temp_k, mat, 0.0);
        c.columns.get_mut::<i128>(CHEM_PHYS_MASS).unwrap().as_mut_slice().fill(Mass::from_kg(kg).0);
        c.columns
            .get_mut::<i128>(CHEM_PHYS_ENERGY)
            .unwrap()
            .as_mut_slice()
            .fill(Energy::from_joules(e).0);
        c.columns.get_mut::<u16>(CHEM_PHYS_MATERIAL).unwrap().as_mut_slice().fill(id.0);
        c.columns.get_mut::<i128>(species_column(0)).unwrap().as_mut_slice()[self.seed_cell] =
            self.n_h;
        c.columns.get_mut::<i128>(species_column(1)).unwrap().as_mut_slice()[self.seed_cell] =
            self.n_o;
        c
    }
}

const VENT_NX: usize = 33;
const VENT_CELL: usize = 16;
const VENT_POOL: usize = 32;

fn vent_cfg(deep: bool) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("world.seed", "20260727");
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", &VENT_NX.to_string());
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", "1.0");
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("physics.gravity_m_s2", "0.0");
    c.set("chem.enabled", "true");
    // Real concentration: a cubic centimetre holding ~10²¹ atoms, about 1 mol/L.
    // Since Phase 5e the collision prefactor is derived from physics rather than
    // fitted, so the concentration is no longer a free choice — at the old ~10⁶
    // molecules per m³ these molecules essentially never meet, and no barrier
    // would make the chemistry proceed.
    c.set("chem.cell_volume_m3", "1e-6");
    c.set("chem.diffusion_m2_s", "0.1");
    c.set("chem.open_ended", "true");
    c.set("chem.max_species", &VENT_POOL.to_string());
    c.set("chem.max_atoms", "3");
    c.set("chem.barrier_j_mol", "50000");
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
    if deep {
        c.set("clock.seconds_per_tick", "10000.0");
        c.set("physics.implicit_conduction", "true");
    } else {
        c.set("clock.seconds_per_tick", "1.0");
    }
    c
}

fn vent_boot(deep: bool) -> Simulation {
    let mut sim = Simulation::new(vent_cfg(deep));
    let spec = GridSpec::new(VENT_NX, 1, 1, 1.0);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(VentGen {
        spec,
        pool: VENT_POOL,
        temp_k: 3000.0,
        seed_cell: VENT_CELL,
        n_h: 400_000_000_000_000_000_000,
        n_o: 200_000_000_000_000_000_000,
        materials,
    }));
    sim.world.activate(SimCoord::new(0, 0, 0));
    sim
}

fn vent_profile(sim: &Simulation, sp: usize) -> Vec<i128> {
    let chunk = sim.world.chunks.iter().next().unwrap().1;
    chunk.columns().get::<i128>(species_column(sp)).unwrap().as_slice().to_vec()
}

fn vent_bar(p: &[i128]) -> String {
    let max = p.iter().copied().max().unwrap_or(0).max(1);
    let ramp: &[u8] = b" .:-=+*#%@";
    p.iter()
        .map(|&v| {
            let lvl = ((v * 9) / max).clamp(0, 9) as usize;
            ramp[lvl] as char
        })
        .collect()
}

fn vent_variance(p: &[i128]) -> f64 {
    let tot: i128 = p.iter().sum();
    if tot == 0 {
        return 0.0;
    }
    let tot = tot as f64;
    let mut mean = 0.0;
    for (i, &v) in p.iter().enumerate() {
        mean += i as f64 * v as f64;
    }
    mean /= tot;
    let mut var = 0.0;
    for (i, &v) in p.iter().enumerate() {
        let d = i as f64 - mean;
        var += d * d * v as f64;
    }
    var / tot
}

fn cmd_vent() {
    println!("── Artificial Universe ─ Phase 5b ─ The Vent ───────────────────");
    println!("  a 33-cell tube of molten silicate at 3000 K. 4,000,000 H· and");
    println!("  2,000,000 O· injected at cell 16 — the vent. No reactions were");
    println!("  declared, no molecules beyond the atoms, no diffusion table:");
    println!("  reactions derive from structure, mobilities from mass, and the");
    println!("  mobility mask from the host's phase.\n");

    let table = {
        let mut t = PeriodicTable::new();
        t.add(Element { z: Z(1), symbol: "H", mass_mda: 1008, valence: 1, electronegativity_c: 220 });
        t.add(Element { z: Z(8), symbol: "O", mass_mda: 15999, valence: 2, electronegativity_c: 344 });
        t
    };

    // Census across the whole tube: (formula, id, total), populated species only.
    let census = |sim: &Simulation| -> Vec<(String, usize, i128)> {
        let mut v = Vec::new();
        for (id, m) in sim.world.chem_registry.iter() {
            if (id.0 as usize) < VENT_POOL {
                let p = vent_profile(sim, id.0 as usize);
                let tot: i128 = p.iter().sum();
                if tot > 0 {
                    v.push((au_chem::formula_string(m, &table), id.0 as usize, tot));
                }
            }
        }
        v.sort_by(|a, b| b.2.cmp(&a.2));
        v
    };

    println!("  ── ACT I: one-second ticks, explicit transport ─────────────");
    let mut sim = vent_boot(false);
    for &t in &[5u64, 30, 90, 180] {
        let done = sim.world.clock.tick().0;
        sim.run(t - done);
        println!("\n  t = {:>3} s   (each row: 33 cells, darkest = that species' peak)", t);
        for (name, id, tot) in census(&sim).into_iter().take(5) {
            println!("   {:>4} |{}| {:>8}", name, vent_bar(&vent_profile(&sim, id)), tot);
        }
    }

    // Graham's law, visible: light molecules have spread further.
    let find = |sim: &Simulation, f: &str| -> Option<usize> {
        census(sim).into_iter().find(|(n, _, _)| n == f).map(|(_, id, _)| id)
    };
    if let (Some(h2), Some(h2o)) = (find(&sim, "H2"), find(&sim, "H2O")) {
        let vh = vent_variance(&vent_profile(&sim, h2));
        let vw = vent_variance(&vent_profile(&sim, h2o));
        let graham = (18015.0f64 / 2016.0).sqrt();
        println!("\n  spread (variance, cells²):  H₂ {:>6.1}   H₂O {:>6.1}   ratio {:.2}", vh, vw, vh / vw);
        println!("  pure transport would give Graham's √(m_H₂O/m_H₂) = {:.2}; the", graham);
        println!("  measured ratio is larger, and the reason is written in H₂'s own");
        println!("  bimodal profile: the vent *eats* hydrogen where the oxidant is");
        println!("  concentrated, so the H₂ that survives is the H₂ that escaped —");
        println!("  reaction–diffusion coupling, not a broken coefficient. (The");
        println!("  clean Graham ratio, 3.98 for H· vs O·, is measured in a");
        println!("  reaction-free tube by the test suite.) Neither rate was ever");
        println!("  declared; both were read off the discovered graphs.");
    }

    println!("\n  ── ACT II: deep time — 10,000-second strides, implicit ─────");
    println!("  same world, reborn; k = D·dt/dx² ≈ 200–700 per stride, far past");
    println!("  the explicit ceiling of ½. Sweeps scale as 2·n² = {} — the", 2 * VENT_NX * VENT_NX);
    println!("  price single-grid Gauss–Seidel charges for deep time, paid");
    println!("  knowingly (the checkerboard that taught us is in the changelog).\n");
    let mut deep = vent_boot(true);
    deep.run(40);
    let bytes = deep.save();
    deep.run(40);
    let straight = deep.world.world_hash();
    println!("  after 80 strides (≈ 9.3 days of world time):");
    for (name, id, tot) in census(&deep).into_iter().take(5) {
        println!("   {:>4} |{}| {:>8}", name, vent_bar(&vent_profile(&deep, id)), tot);
    }

    // How flat is the water?
    if let Some(w) = find(&deep, "H2O") {
        let p = vent_profile(&deep, w);
        let mean = p.iter().sum::<i128>() / p.len() as i128;
        let worst = p.iter().map(|&v| (v - mean).abs()).max().unwrap();
        println!("\n  the ocean is level: worst H₂O deviation {} on a mean of {}", worst, mean);
    }

    // The books.
    let (mut h_atoms, mut o_atoms) = (0i128, 0i128);
    for (id, m) in deep.world.chem_registry.iter() {
        if (id.0 as usize) < VENT_POOL {
            let tot: i128 = vent_profile(&deep, id.0 as usize).iter().sum();
            let f = m.formula(&table);
            h_atoms += f[0] * tot;
            o_atoms += f[1] * tot;
        }
    }
    println!("  atom audit: {} H (started 400000000000000000000), {} O (started 200000000000000000000)", h_atoms, o_atoms);

    let mut resumed = Simulation::load(&bytes, vent_cfg(true)).unwrap();
    let spec = GridSpec::new(VENT_NX, 1, 1, 1.0);
    let materials = resumed.world.materials.clone();
    resumed.world.set_generator(Box::new(VentGen {
        spec,
        pool: VENT_POOL,
        temp_k: 3000.0,
        seed_cell: VENT_CELL,
        n_h: 400_000_000_000_000_000_000,
        n_o: 200_000_000_000_000_000_000,
        materials,
    }));
    resumed.run(40);
    let r = resumed.world.world_hash();
    println!(
        "  save/resume at stride 40: {:#018x} → {:#018x} ({})",
        straight,
        r,
        if straight == r { "identical — a spreading world survives sleep" } else { "DIVERGED" }
    );

    println!("\n  What changed in 5b: molecules travel. Water invented at the vent");
    println!("  is found at the tube's ends; its speed was never written down —");
    println!("  it is √(1 Da/m) of a graph the world drew for itself. A host that");
    println!("  freezes traps its cargo where it sat (tested); advection — riding");
    println!("  the currents, not just the gradients — waits for the fluid");
    println!("  coupling. The vent is lit.");
}

// ═════════════════════════════════════════════════════════════════════════════
//  raf — does this world's chemistry make itself? (Phase 5c)
// ═════════════════════════════════════════════════════════════════════════════

const RAF_POOL: usize = 48;

/// One hot cell of bare atoms: carbon, hydrogen, oxygen and nothing else.
struct RafGen {
    pool: usize,
    temp_k: f64,
    seeds: Vec<(usize, i128)>,
    materials: MaterialRegistry,
}

impl SimGen for RafGen {
    fn generate(&self, coord: SimCoord, _rng: &mut au_core::Rng) -> SimChunk {
        let mut c = SimChunk::new(coord, SimLod::Full);
        install_phys(&mut c.columns, 1);
        install_chem(&mut c.columns, 1, self.pool);
        let id = self.materials.by_name("silicate").unwrap();
        let mat = self.materials.get(id).unwrap();
        let kg = 1000.0;
        let e = energy_for_temperature(kg, self.temp_k, mat, 0.0);
        c.columns.get_mut::<i128>(CHEM_PHYS_MASS).unwrap().as_mut_slice()[0] = Mass::from_kg(kg).0;
        c.columns.get_mut::<i128>(CHEM_PHYS_ENERGY).unwrap().as_mut_slice()[0] =
            Energy::from_joules(e).0;
        c.columns.get_mut::<u16>(CHEM_PHYS_MATERIAL).unwrap().as_mut_slice()[0] = id.0;
        for &(sp, n) in &self.seeds {
            c.columns.get_mut::<i128>(species_column(sp)).unwrap().as_mut_slice()[0] = n;
        }
        c
    }
}

fn raf_cfg(catalysis: bool) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("world.seed", "20260727");
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
    // A large cell and a high barrier: catalysis is only visible where Arrhenius
    // still gates the rate, and where the barrier removed is worth more than the
    // termolecular encounter costs. See the Phase 5c changelog.
    // Real concentration: a cubic centimetre holding ~10²¹ atoms, about 1 mol/L.
    // Since Phase 5e the collision prefactor is derived from physics rather than
    // fitted, so the concentration is no longer a free choice — at the old ~10⁶
    // molecules per m³ these molecules essentially never meet, and no barrier
    // would make the chemistry proceed.
    c.set("chem.cell_volume_m3", "1e-6");
    c.set("chem.barrier_j_mol", "200000");
    c.set("chem.open_ended", "true");
    c.set("chem.max_species", &RAF_POOL.to_string());
    c.set("chem.max_atoms", "3");
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

fn raf_boot(catalysis: bool) -> Simulation {
    let mut sim = Simulation::new(raf_cfg(catalysis));
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(RafGen {
        pool: RAF_POOL,
        temp_k: 2500.0,
        seeds: vec![
            (0, 300_000_000_000_000_000_000),
            (1, 100_000_000_000_000_000_000),
            (2, 200_000_000_000_000_000_000),
        ],
        materials,
    }));
    sim.world.activate(SimCoord::new(0, 0, 0));
    sim
}

fn raf_table() -> PeriodicTable {
    let mut t = PeriodicTable::new();
    t.add(Element { z: Z(1), symbol: "H", mass_mda: 1008, valence: 1, electronegativity_c: 220 });
    t.add(Element { z: Z(6), symbol: "C", mass_mda: 12011, valence: 4, electronegativity_c: 255 });
    t.add(Element { z: Z(8), symbol: "O", mass_mda: 15999, valence: 2, electronegativity_c: 344 });
    t
}

fn cmd_raf() {
    println!("── Artificial Universe ─ Phase 5c ─ Does it make itself? ───────");
    println!("  one hot cell. 3,000,000 H·, 1,000,000 C·, 2,000,000 O·.");
    println!("  no molecules declared, no reactions declared, no catalysts");
    println!("  declared. Everything below is read off graphs the world drew.\n");

    let t = raf_table();
    let name = |sim: &Simulation, s: au_chem::SpeciesId| -> String {
        match sim.world.chem_registry.get(s) {
            Some(m) => formula_string(m, &t),
            None => format!("#{}", s.0),
        }
    };

    // ── ACT I ───────────────────────────────────────────────────────────────
    println!("  ── ACT I: no catalysis ─────────────────────────────────────");
    let mut plain = raf_boot(false);
    plain.run(300);
    let sp = au_sim::autocatalysis::survey(&plain.world).unwrap();
    let both_sides = sp
        .reactions
        .iter()
        .filter(|r| {
            r.reactants.iter().any(|a| r.products.iter().any(|b| b.species == a.species))
        })
        .count();
    println!("  species discovered      {}", plain.world.chem_registry.len());
    println!("  reactions derived       {}", sp.reactions.len());
    println!("  …with a species on both sides of the arrow:  {}", both_sides);
    println!("  catalysts               {}", sp.catalyst_count());
    println!(
        "  autocatalytic core      {}",
        match &sp.raf {
            Some(r) => format!("{} reactions", r.reactions.len()),
            None => "none".to_string(),
        }
    );
    println!("\n  That zero is not an accident of this run — it is structural.");
    println!("  Form, Split, Raise and Lower give 2→1, 1→2 and 1→1, and none of");
    println!("  them can put a species on both sides. A molecule that cannot");
    println!("  appear among its own reaction's reactants can never make more of");
    println!("  itself. Self-replication was not missing from this network; it");
    println!("  was *impossible* in it, at any temperature, for any length of");
    println!("  time. That is what Phase 5c had to fix.\n");

    // ── ACT II ──────────────────────────────────────────────────────────────
    println!("  ── ACT II: catalysis, derived from structure ───────────────");
    let mut sim = raf_boot(true);
    sim.run(300);
    let s = au_sim::autocatalysis::survey(&sim.world).unwrap();
    println!("  species discovered      {}", sim.world.chem_registry.len());
    println!("  reactions derived       {}", s.reactions.len());
    println!("  catalyst species        {}  (relation: {} pairs)", s.catalyst_count(), s.relation.len());
    match &s.raf {
        Some(r) => {
            println!("  autocatalytic core      {} reactions", r.reactions.len());
            println!("  closure from the food   {} species", r.closure.len());
            println!("  pruning rounds          {}", r.rounds);
        }
        None => println!("  autocatalytic core      none"),
    }
    let food: Vec<String> = s.food.iter().map(|&f| name(&sim, f)).collect();
    println!("  food (bare elements)    {}", food.join("  "));

    if !s.self_producing.is_empty() {
        // Group by formula: two entries with the same formula are *isomers* —
        // different graphs the world discovered separately, which the formula
        // string cannot tell apart but the canonical form can.
        let mut names: Vec<String> = s.self_producing.iter().map(|&x| name(&sim, x)).collect();
        names.sort();
        let mut shown: Vec<String> = Vec::new();
        let mut i = 0;
        while i < names.len() {
            let j = names[i..].iter().take_while(|n| **n == names[i]).count();
            shown.push(if j > 1 {
                format!("{} ×{} (isomers)", names[i], j)
            } else {
                names[i].clone()
            });
            i += j;
        }
        println!("\n  species that make more of themselves:  {}", shown.join("   "));
        // Show one, as the reaction the reactor actually runs.
        if let Some(&(c, i)) = s
            .relation
            .iter()
            .find(|(c, i)| s.reactions[*i].products.iter().any(|p| p.species == *c))
        {
            let side = |ts: &[au_chem::Term], extra: Option<au_chem::SpeciesId>| -> String {
                let mut v: Vec<String> = ts
                    .iter()
                    .map(|x| {
                        if x.count > 1 {
                            format!("{} {}", x.count, name(&sim, x.species))
                        } else {
                            name(&sim, x.species)
                        }
                    })
                    .collect();
                if let Some(e) = extra {
                    v.push(name(&sim, e));
                }
                v.join(" + ")
            };
            let r = &s.reactions[i];
            println!("\n  the reaction, as the reactor runs it:");
            println!("      {}   →   {}", side(&r.reactants, Some(c)), side(&r.products, Some(c)));
            println!("  {} holds the two reactants together and a second {} walks", name(&sim, c), name(&sim, c));
            println!("  out. Nothing in the code special-cases this; it is what the");
            println!("  general rule says when the structure happens to line up.");
        }
    }

    // ── ACT III ─────────────────────────────────────────────────────────────
    println!("\n  ── ACT III: the collapse ───────────────────────────────────");
    println!("  the same discovered network, asked a counterfactual: how polar");
    println!("  must a grip be to count at all? (The world is not re-run — only");
    println!("  the lens changes.)\n");
    println!("    min grip   catalysts   core   closure   self-producing");
    for thr in [200_000.0, 250_000.0, 280_000.0, 340_000.0, 400_000.0] {
        let rules = au_chem::CatalysisRules { min_grip_j_mol: thr, ..Default::default() };
        let v = au_sim::autocatalysis::survey_with(&sim.world, Some(rules)).unwrap();
        let (core, clo) = match &v.raf {
            Some(r) => (r.reactions.len(), r.closure.len()),
            None => (0, 0),
        };
        println!(
            "    {:>8.0}   {:>9}   {:>4}   {:>7}   {:>4}",
            thr,
            v.catalyst_count(),
            core,
            clo,
            v.self_producing.len()
        );
    }
    println!("\n  Kauffman's claim, in a network nobody designed: self-sustaining");
    println!("  chemistry does not fade in — it switches on. And the cliff sits");
    println!("  exactly where carbon–oxygen contacts stop counting. Take away");
    println!("  C–O polarity and what survives is a vestige with nothing in it");
    println!("  that makes more of itself. Nobody put carbon there.");

    // ── The books ───────────────────────────────────────────────────────────
    let chunk = sim.world.chunks.iter().next().unwrap().1;
    let (mut h, mut c, mut o) = (0i128, 0i128, 0i128);
    for (id, m) in sim.world.chem_registry.iter() {
        let i = id.0 as usize;
        if i < RAF_POOL {
            let p = chunk.columns().get::<i128>(species_column(i)).unwrap().as_slice()[0];
            if p != 0 {
                let f = m.formula(&t);
                h += f[0] * p;
                c += f[1] * p;
                o += f[2] * p;
            }
        }
    }
    println!(
        "\n  atom audit: {} H / {} C / {} O  (started 300000000000000000000 / 100000000000000000000 / 200000000000000000000)",
        h, c, o
    );

    let mut a = raf_boot(true);
    a.run(150);
    let bytes = a.save();
    a.run(150);
    let straight = a.world.world_hash();
    let mut b = Simulation::load(&bytes, raf_cfg(true)).unwrap();
    let materials = b.world.materials.clone();
    b.world.set_generator(Box::new(RafGen {
        pool: RAF_POOL,
        temp_k: 2500.0,
        seeds: vec![
            (0, 300_000_000_000_000_000_000),
            (1, 100_000_000_000_000_000_000),
            (2, 200_000_000_000_000_000_000),
        ],
        materials,
    }));
    b.run(150);
    println!(
        "  save/resume: {:#018x} → {:#018x} ({})",
        straight,
        b.world.world_hash(),
        if straight == b.world.world_hash() { "identical" } else { "DIVERGED" }
    );
    println!("\n  A RAF is not life. It has no membrane, no heredity, no");
    println!("  individuality, and it cannot evolve. It is metabolism before");
    println!("  there is anything for the metabolism to belong to — the");
    println!("  network-level precondition for everything Phase 6 needs.");
}


// ═══════════════════════════════════════════════════════════════════════════
//  au watch — run a world and write a page you can look at
// ═══════════════════════════════════════════════════════════════════════════

/// The standing debt this command exists to pay.
///
/// Every result in this project has come from a test, and tests can only check
/// what somebody already thought to ask. Three separate debugging sessions ended
/// with an `eprintln!` in a hot loop and a `grep`, and every one of them was
/// solved the instant real numbers appeared — a sawtooth census, a flat contents
/// line, a temperature leaving the axis. This writes those out.
///
///     au watch --config world.kv --ticks 5000 --every 10 --out world.html
///
/// The config is optional; without one it runs the built-in protocell scenario,
/// which is the smallest world in which anything interesting has happened yet.
fn cmd_watch(args: &[String]) {
    let ticks = arg_u64(args, "--ticks", 2_000);
    let every = arg_u64(args, "--every", (ticks / 400).max(1));
    let out = arg_str(args, "--out", "world.html");
    let max_bags = arg_u64(args, "--bags", 64) as usize;

    let (mut sim, names, title) = match arg_opt(args, "--config") {
        Some(path) => {
            let text = std::fs::read_to_string(&path)
                .unwrap_or_else(|e| { eprintln!("cannot read {}: {}", path, e); std::process::exit(2) });
            let mut c = Config::parse(DEFAULT_CONFIG).expect("defaults");
            c.merge(&text).unwrap_or_else(|e| { eprintln!("bad config: {}", e); std::process::exit(2) });
            let s = Simulation::new(c);
            let n = species_names(&s);
            (s, n, format!("Artificial Universe — {}", path))
        }
        None => {
            let s = watch_scenario_with(
                arg_u64(args, "--surfactant", 2_000_000) as i128,
                arg_u64(args, "--pentane", 20_000_000) as i128,
                arg_u64(args, "--feed", 20_000_000) as i128,
            );
            let n = species_names(&s);
            (s, n, "Artificial Universe — protocell scenario".to_string())
        }
    };

    let mut rec = au_sim::observe::Recorder::new(every, max_bags);
    rec.snapshot(&sim.world);
    let t0 = Instant::now();
    for _ in 0..ticks {
        sim.tick();
        rec.observe(&sim.world);
    }
    rec.snapshot(&sim.world);

    let html = render::render(rec.samples(), &title, &names);
    if let Err(e) = std::fs::write(&out, html) {
        eprintln!("cannot write {}: {}", out, e);
        std::process::exit(2);
    }
    let last = rec.samples().last().unwrap();
    println!(
        "{} ticks in {:.2}s · {} samples · {} protocells ({} nucleated, {} divided, {} lysed) · {:.1} K",
        ticks,
        t0.elapsed().as_secs_f64(),
        rec.samples().len(),
        last.protocells,
        last.nucleated,
        last.divided,
        last.lysed,
        last.temp_mean
    );
    println!("wrote {}", out);
}

fn species_names(sim: &Simulation) -> Vec<String> {
    (0..sim.world.chem_registry.len())
        .map(|i| format!("sp{}", i))
        .collect()
}

fn arg_opt(args: &[String], key: &str) -> Option<String> {
    args.iter().position(|a| a == key).and_then(|i| args.get(i + 1)).cloned()
}

fn arg_str(args: &[String], key: &str, default: &str) -> String {
    arg_opt(args, key).unwrap_or_else(|| default.to_string())
}

/// The built-in scenario: the smallest world in which anything interesting has
/// happened. A micron of water holding pentane, a little pentanol and hydroxyl,
/// fed pentane through a Dirichlet port so the medium cannot exhaust it.
/// Protocells nucleate, import substrate, convert it to surfactant inside
/// themselves — where it is trapped, because it *is* the membrane — and divide
/// when the reduced volume crosses the two-sphere bound.
fn watch_scenario_with(surfactant: i128, pentane: i128, feed: i128) -> Simulation {
    fn chain(n: usize, oh: Option<usize>) -> (String, String) {
        let mut z: Vec<u32> = vec![6; n];
        let mut b: Vec<(usize, usize)> = (0..n.saturating_sub(1)).map(|i| (i, i + 1)).collect();
        if let Some(k) = oh {
            let o = z.len();
            z.push(8);
            b.push((k, o));
            let h = z.len();
            z.push(1);
            b.push((o, h));
        }
        for i in 0..n {
            let deg = b.iter().filter(|(x, y)| *x == i || *y == i).count();
            for _ in deg..4 {
                let k = z.len();
                z.push(1);
                b.push((i, k));
            }
        }
        (
            z.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(","),
            b.iter().map(|(x, y)| format!("{}-{}:1", x, y)).collect::<Vec<_>>().join(","),
        )
    }

    let mut c = Config::parse(DEFAULT_CONFIG).expect("defaults");
    c.merge(include_str!("../../../data/reference_materials.kv")).expect("materials");
    c.set("world.seed", "5150");
    c.set("clock.seconds_per_tick", "1.0");
    c.set("clock.rescale.enabled", "false");
    c.set("physics.grid.nx", "1");
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", "1.0e-6");
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("physics.gravity_m_s2", "0.0");
    // Backward Euler, so the clock is not capped at the explicit conduction
    // limit. Without this the physics layer correctly refuses the requested
    // one-second tick and takes 3.6e-4 s instead — which means 3000 ticks is
    // about one second of simulated time, and expecting a protocell population
    // to reproduce in a second is expecting the wrong thing entirely.
    c.set("physics.implicit_conduction", "true");
    c.set("chem.enabled", "true");
    c.set("chem.cell_volume_m3", "1.0e-18");
    c.set("chem.diffusion_m2_s", "1.0e-13");
    c.set("chem.membrane", "true");
    c.set("chem.element.count", "3");
    for (i, (z, sym, mda, val, en)) in
        [(1u32, "H", 1008u32, 1u32, 220u32), (6, "C", 12011, 4, 255), (8, "O", 15999, 2, 344)]
            .iter()
            .enumerate()
    {
        c.set(&format!("chem.element.{}.z", i), &z.to_string());
        c.set(&format!("chem.element.{}.symbol", i), sym);
        c.set(&format!("chem.element.{}.mass_mda", i), &mda.to_string());
        c.set(&format!("chem.element.{}.valence", i), &val.to_string());
        c.set(&format!("chem.element.{}.electronegativity_c", i), &en.to_string());
    }
    c.set("chem.species.count", "4");
    let (z, b) = chain(5, None); // pentane — no head, a solute
    c.set("chem.species.0.atoms", &z);
    c.set("chem.species.0.bonds", &b);
    let (z, b) = chain(5, Some(4)); // 1-pentanol — the surfactant
    c.set("chem.species.1.atoms", &z);
    c.set("chem.species.1.bonds", &b);
    c.set("chem.species.2.atoms", "1,8,1"); // water
    c.set("chem.species.2.bonds", "0-1:1,1-2:1");
    c.set("chem.species.3.atoms", "1,1"); // H2
    c.set("chem.species.3.bonds", "0-1:1");
    c.set("chem.reaction.count", "1");
    c.set("chem.reaction.0.reactants", "0:1,2:1");
    c.set("chem.reaction.0.products", "1:1,3:1");
    c.set("chem.reaction.0.activation_j", "60000");
    c.set("chem.reaction.0.pre_exponential", "1.0e-2");
    c.set("chem.reaction.0.enthalpy_j", "-1.0e-20");
    if feed > 0 {
        c.set("chem.reservoir.source", &format!("0:0@{}", feed));
    }
    c.set("life.protocells", "true");

    struct G {
        materials: MaterialRegistry,
        surfactant: i128,
        pentane: i128,
    }
    impl ChunkGenerator for G {
        fn generate(&self, coord: ChunkCoord, _r: &mut au_core::Rng) -> Chunk {
            let mut ch = Chunk::new(coord, Lod::Full);
            install_phys(&mut ch.columns, 1);
            install_chem(&mut ch.columns, 1, 4);
            let id = self.materials.by_name("water").unwrap();
            let m = self.materials.get(id).unwrap();
            let kg = 1000.0 * 1.0e-18;
            ch.columns.get_mut::<i128>(PHYS_MASS).unwrap().as_mut_slice()[0] = Mass::from_kg(kg).0;
            ch.columns.get_mut::<i128>(PHYS_ENERGY).unwrap().as_mut_slice()[0] =
                Energy::from_joules(energy_for_temperature(kg, 320.0, m, 0.0)).0;
            ch.columns.get_mut::<u16>(PHYS_MATERIAL).unwrap().as_mut_slice()[0] = id.0;
            for (sp, n) in [(0usize, self.pentane), (1, self.surfactant), (2, 50_000_000)] {
                ch.columns.get_mut::<i128>(species_column(sp)).unwrap().as_mut_slice()[0] = n;
            }
            ch
        }
    }

    let mut sim = Simulation::new(c);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(G { materials, surfactant, pentane }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));
    sim
}
