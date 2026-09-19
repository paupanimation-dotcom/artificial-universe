//! Species transport in the world — Phase 5b's integration.
//!
//! au-physics proved the diffusion operator (conservation, Einstein's variance,
//! the implicit twin, the donor's-pocket limiter). These tests prove the whole
//! tower: a species' diffusion rate is derived from the mass of its *discovered
//! molecular graph*, the mobility mask is derived from the host's phase by the
//! same `derive()` everything trusts, and chemistry + transport together keep
//! atoms exact across cells, keep determinism, and survive resume.
//!
//! The anchor is **Graham's law** (1846): diffusion rate ∝ 1/√m. Hydrogen and
//! oxygen atoms differ in mass by ×15.87, so their spreads must differ by
//! ×√15.87 ≈ 3.98. To be clear about what is being tested: the √m scaling is
//! *built into* the coefficient model — the law here is input, not output. What
//! the test validates is the tower's plumbing: that the mass is read correctly
//! off the graph, the per-species rates thread through the scheduler, the
//! transport delivers them unmangled, and the measured variances land on the
//! analytic ratio end-to-end. (Ra_c played the same role for the fluid solver.)

use au_core::hash::WorldHash;
use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, Lod};
use au_physics::{energy_for_temperature, Energy, GridSpec, Mass};
use au_sim::chem_columns::{install as install_chem, species_column};
use au_sim::physics_columns::{install as install_phys, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};
use au_sim::{Config, Simulation, DEFAULT_CONFIG};

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");

const NX: usize = 41;
const CENTER: usize = 20;

/// A 1-D tube of host material, with species seeded in the centre cell only.
struct TubeGen {
    spec: GridSpec,
    pool: usize,
    temp_k: f64,
    seed_cell: usize,
    seeds: Vec<(usize, i128)>, // (species column, count) — placed in seed_cell
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
        for &(sp, count) in &self.seeds {
            c.columns.get_mut::<i128>(species_column(sp)).unwrap().as_mut_slice()
                [self.seed_cell] = count;
        }
        c
    }
}

/// A tube world with hydrogen and oxygen atoms declared and NO reactions —
/// pure transport, so Graham's ratio is unpolluted by chemistry.
fn tube_config(seed: u64, temp_note: &str) -> Config {
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
    c.set("chem.diffusion_m2_s", "0.1");
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
    let _ = temp_note;
    c
}

fn tube_sim(seed: u64, temp_k: f64) -> Simulation {
    let mut sim = Simulation::new(tube_config(seed, ""));
    let spec = GridSpec::new(NX, 1, 1, 1.0);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(TubeGen {
        spec,
        pool: 2,
        temp_k,
        seed_cell: CENTER,
        seeds: vec![(0, 1_000_000_000), (1, 1_000_000_000)],
        materials,
    }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));
    sim
}

fn profile(sim: &Simulation, sp: usize) -> Vec<i128> {
    let chunk = sim.world.chunks.iter().next().unwrap().1;
    chunk.columns().get::<i128>(species_column(sp)).unwrap().as_slice().to_vec()
}

fn variance(p: &[i128]) -> f64 {
    let tot: i128 = p.iter().sum();
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

/// **Graham's law, end to end.** Equal releases of H· and O· in a liquid tube;
/// after 400 ticks, the variance ratio must be √(m_O/m_H) = √(15999/1008)
/// ≈ 3.984 — the per-species rates were read off the molecular graphs and
/// delivered through the scheduler intact.
/// (Measured at 150 ticks: hydrogen's σ is then ~5.5 cells against walls
/// 20 cells out — free-space Gaussians. Run longer and the sealed ends start
/// reflecting the fast species' tail, saturating its variance toward the
/// uniform limit and dragging the ratio down: a finite-tube truth, not a bug,
/// and the first version of this test learned it at 400 ticks.)
#[test]
fn hydrogen_outruns_oxygen_by_grahams_ratio() {
    let mut sim = tube_sim(11, 1500.0); // silicate melts at 1400 K: a magma ocean
    sim.run(150);

    let h = profile(&sim, 0);
    let o = profile(&sim, 1);
    assert_eq!(h.iter().sum::<i128>(), 1_000_000_000, "hydrogen leaked");
    assert_eq!(o.iter().sum::<i128>(), 1_000_000_000, "oxygen leaked");

    let ratio = variance(&h) / variance(&o);
    let graham = (15999.0f64 / 1008.0).sqrt();
    assert!(
        (ratio / graham - 1.0).abs() < 0.04,
        "spread ratio {:.3} vs Graham's {:.3} — the mass-derived rates were mangled in transit",
        ratio,
        graham
    );
    // And hydrogen has genuinely travelled: the tube edges have seen it.
    assert!(h[2] > 0, "hydrogen never reached the tube's reaches");
}

/// A cold tube is a solid tube, and a solid host is a wall everywhere: the
/// species sit exactly where they were born, to the last count. This is the
/// mask being derived from phase — the same freeze that stops an ocean stops
/// its chemistry.
#[test]
fn a_frozen_host_traps_its_chemistry_in_place() {
    let mut sim = tube_sim(12, 300.0); // solid silicate
    sim.run(400);
    let h = profile(&sim, 0);
    let o = profile(&sim, 1);
    for i in 0..NX {
        let expect = if i == CENTER { 1_000_000_000 } else { 0 };
        assert_eq!(h[i], expect, "hydrogen moved through solid rock at cell {}", i);
        assert_eq!(o[i], expect, "oxygen moved through solid rock at cell {}", i);
    }
}

/// Deep time for matter: the same tube with `physics.implicit_conduction = true`
/// and ticks ten thousand seconds long — k = D·dt/dx² ≈ 1000, absurdly past the
/// explicit ceiling. It runs, it conserves to the last count, and two runs of
/// the same seed agree to the bit.
///
/// At this k the donor's-pocket limiter governs the pace: each stride can move
/// at most half a cell's holdings per face, so mixing proceeds at roughly an
/// explicit k of 0.5 per stride regardless of how huge k is — stability and
/// non-negativity bought at the price of no teleportation. 200 strides is
/// comfortably past the tube's saturation on that clock.
#[test]
fn deep_time_diffusion_is_stable_conserved_and_deterministic() {
    let run = || {
        let mut cfg = tube_config(13, "");
        cfg.set("clock.seconds_per_tick", "10000.0");
        cfg.set("physics.implicit_conduction", "true");
        let mut sim = Simulation::new(cfg);
        let spec = GridSpec::new(NX, 1, 1, 1.0);
        let materials = sim.world.materials.clone();
        sim.world.set_generator(Box::new(TubeGen {
            spec,
            pool: 2,
            temp_k: 1500.0,
            seed_cell: CENTER,
            seeds: vec![(0, 1_000_000_000), (1, 1_000_000_000)],
            materials,
        }));
        sim.world.activate(ChunkCoord::new(0, 0, 0));
        sim.run(200);
        let h = profile(&sim, 0);
        let o = profile(&sim, 1);
        assert_eq!(h.iter().sum::<i128>(), 1_000_000_000, "deep time leaked hydrogen");
        assert_eq!(o.iter().sum::<i128>(), 1_000_000_000, "deep time leaked oxygen");
        assert!(h.iter().all(|&v| v >= 0) && o.iter().all(|&v| v >= 0), "negative population");
        // At k ≈ 1000 the pocket limiter sets the pace; 200 strides is past
        // saturation on that clock.
        let mean = 1_000_000_000 / NX as i128;
        let worst = h.iter().map(|&v| (v - mean).abs()).max().unwrap();
        assert!(
            worst < mean / 2,
            "deep-time hydrogen failed to mix: worst deviation {} vs mean {} (profile {:?})",
            worst,
            mean,
            h
        );
        sim.world.world_hash()
    };
    assert_eq!(run(), run(), "deep-time diffusion lost determinism");
}

/// Genesis meets transport: an open chemistry seeded only at the centre
/// discovers water *and* the water travels — found in cells its ingredients
/// were never placed in — while the atom books stay exact across the tube.
#[test]
fn discovered_water_spreads_beyond_its_birthplace() {
    let mut cfg = tube_config(14, "");
    cfg.set("physics.grid.nx", "9");
    cfg.set("chem.open_ended", "true");
    cfg.set("chem.max_species", "32");
    cfg.set("chem.max_atoms", "3");
    cfg.set("chem.barrier_j_mol", "50000");
    // A real concentration. At the tube's default of ~10⁶ molecules per m³ —
    // about 10⁻¹⁷ mol/L — no barrier whatsoever would let this chemistry
    // proceed: even collision-limited, the encounters are too rare. That world
    // was never chemically viable; it only looked viable while the prefactor was
    // ~10¹⁹ times the gas-kinetic value. Here the cell is a cubic centimetre
    // holding ~10²¹ atoms, and the chemistry runs to completion quickly while
    // transport spreads the product slowly — which is the honest ordering of
    // these two processes at any real concentration.
    cfg.set("chem.cell_volume_m3", "1e-6");
    let mut sim = Simulation::new(cfg);
    let spec = GridSpec::new(9, 1, 1, 1.0);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(TubeGen {
        spec,
        pool: 32,
        temp_k: 3000.0, // hot enough to anneal past the hydroxyl trap
        seed_cell: 4,
        seeds: vec![(0, 400_000_000_000_000_000_000), (1, 200_000_000_000_000_000_000)],
        materials,
    }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));
    sim.run(300);

    let water = au_chem::Molecule::new(
        vec![
            au_chem::Atom { z: au_chem::Z(1) },
            au_chem::Atom { z: au_chem::Z(8) },
            au_chem::Atom { z: au_chem::Z(1) },
        ],
        vec![
            au_chem::Bond::new(0, 1, au_chem::BondOrder::Single),
            au_chem::Bond::new(1, 2, au_chem::BondOrder::Single),
        ],
    );
    let w = sim.world.chem_registry.lookup(&water).expect("water should be discovered");
    let wp = profile(&sim, w.0 as usize);
    assert!(wp[4] > 0, "no water at the birthplace");
    assert!(
        wp.iter().enumerate().any(|(i, &v)| i != 4 && v > 0),
        "water never left the cell it was invented in"
    );

    // The atom books, across every cell and every discovered species.
    let table = {
        let mut t = au_chem::PeriodicTable::new();
        t.add(au_chem::Element { z: au_chem::Z(1), symbol: "H", mass_mda: 1008, valence: 1, electronegativity_c: 220 });
        t.add(au_chem::Element { z: au_chem::Z(8), symbol: "O", mass_mda: 15999, valence: 2, electronegativity_c: 344 });
        t
    };
    let (mut h, mut o) = (0i128, 0i128);
    for (id, m) in sim.world.chem_registry.iter() {
        if (id.0 as usize) < 32 {
            let tot: i128 = profile(&sim, id.0 as usize).iter().sum();
            let f = m.formula(&table);
            h += f[0] * tot;
            o += f[1] * tot;
        }
    }
    assert_eq!(h, 400_000_000_000_000_000_000, "hydrogen atoms leaked across the tube");
    assert_eq!(o, 200_000_000_000_000_000_000, "oxygen atoms leaked across the tube");
}

/// A world saved mid-spread and resumed equals one that never stopped.
#[test]
fn spreading_survives_save_and_resume() {
    let mut sim = tube_sim(15, 1500.0);
    sim.run(100);
    let bytes = sim.save();
    sim.run(100);
    let straight = sim.world.world_hash();

    let mut resumed = Simulation::load(&bytes, tube_config(15, "")).expect("resume");
    let spec = GridSpec::new(NX, 1, 1, 1.0);
    let materials = resumed.world.materials.clone();
    resumed.world.set_generator(Box::new(TubeGen {
        spec,
        pool: 2,
        temp_k: 1500.0,
        seed_cell: CENTER,
        seeds: vec![(0, 1_000_000_000), (1, 1_000_000_000)],
        materials,
    }));
    resumed.run(100);
    assert_eq!(straight, resumed.world.world_hash(), "a resumed spread diverged");
}
