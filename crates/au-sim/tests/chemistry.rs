//! Chemistry, wired into the real simulation loop — Phase 3b.
//!
//! The au-chem crate already proved, in isolation, that the reactor conserves
//! atoms and energy and reproduces van 't Hoff. What these tests prove is the
//! *integration*: that chemistry runs inside the scheduler, on the active set,
//! feeding the **same** thermal field physics conducts and radiates — and that
//! doing so leaves determinism and save/resume exactly as intact as they were
//! before chemistry existed.
//!
//! The claim that matters here is not "reactions happen" (au-chem showed that) but
//! "reactions happen *in the world*, changing the world's own energy field, and
//! the world remains a deterministic, resumable function of its seed and config".

use au_core::hash::WorldHash;
use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, Lod};
use au_physics::{energy_for_temperature, Energy, GridSpec, Mass, MaterialRegistry};
use au_sim::chem_columns::{install as install_chem, species_column};
use au_sim::physics_columns::{install as install_phys, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};
use au_sim::{Config, Simulation, DEFAULT_CONFIG};

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");

/// A single-cell parcel of reacting gas. An initial condition, like the seed —
/// not content. It installs physics columns (mass, energy, material) and chemistry
/// columns (species populations), then seeds a stoichiometric mix of the two
/// reactants with zero product, at a chosen starting temperature.
struct GasParcelGen {
    spec: GridSpec,
    temp_k: f64,
    /// molecules of each reactant (species 0 and species 1), per cell
    n_reactant: i128,
    n_species: usize,
    rho: f64,
    materials: MaterialRegistry,
}

impl ChunkGenerator for GasParcelGen {
    fn generate(&self, coord: ChunkCoord, _rng: &mut au_core::Rng) -> Chunk {
        let mut c = Chunk::new(coord, Lod::Full);
        let n = self.spec.cells();
        install_phys(&mut c.columns, n);
        install_chem(&mut c.columns, n, self.n_species);

        // A gas material to carry heat capacity. We reuse the fixture's "air" if
        // present, else silicate — the point is only that the cell has a defined
        // specific heat so temperature is derivable from energy.
        let id = self
            .materials
            .by_name("air")
            .or_else(|| self.materials.by_name("silicate"))
            .expect("fixture must define a material");
        let mat = self.materials.get(id).unwrap();
        let kg = self.rho * self.spec.cell_volume();
        let e = energy_for_temperature(kg, self.temp_k, mat, 0.0);

        c.columns.get_mut::<i128>(PHYS_MASS).unwrap().as_mut_slice().fill(Mass::from_kg(kg).0);
        c.columns.get_mut::<i128>(PHYS_ENERGY).unwrap().as_mut_slice().fill(Energy::from_joules(e).0);
        c.columns.get_mut::<u16>(PHYS_MATERIAL).unwrap().as_mut_slice().fill(id.0);

        // Seed reactants (species 0, 1), leave product (species 2) at zero.
        c.columns.get_mut::<i128>(species_column(0)).unwrap().as_mut_slice().fill(self.n_reactant);
        c.columns.get_mut::<i128>(species_column(1)).unwrap().as_mut_slice().fill(self.n_reactant);
        c
    }
}

/// A config declaring the H₂ + O₂ ⇌ 2 H–O chemistry the au-chem demo uses, plus a
/// single-cell physics grid. The reaction is exothermic and given a low barrier so
/// it fires readily at the test temperature — the point is to see heat enter the
/// shared field, not to study kinetics (au-chem did that).
fn reacting_config(seed: u64) -> Config {
    let mut c = Config::parse(DEFAULT_CONFIG).unwrap();
    c.merge(REFERENCE_MATERIALS).unwrap();
    c.set("world.seed", &seed.to_string());
    c.set("clock.seconds_per_tick", "1.0");
    c.set("clock.rescale.enabled", "false");

    // One cell. Chemistry is per-cell; a 1×1×1 grid isolates the coupling.
    c.set("physics.grid.nx", "1");
    c.set("physics.grid.ny", "1");
    c.set("physics.grid.nz", "1");
    c.set("physics.grid.cell_m", "1.0");
    // No boundary heat exchange: we want to see ONLY chemistry move the energy, so
    // physics conducts/radiates nothing in or out.
    c.set("physics.bottom_flux_w_m2", "0.0");
    c.set("physics.top_radiates", "false");
    c.set("physics.gravity_m_s2", "0.0");

    // ── The declared chemistry vocabulary. ──
    c.set("chem.enabled", "true");
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

    c.set("chem.species.count", "3");
    c.set("chem.species.0.atoms", "1,1"); // H2
    c.set("chem.species.0.bonds", "0-1:1");
    c.set("chem.species.1.atoms", "8,8"); // O2
    c.set("chem.species.1.bonds", "0-1:1");
    c.set("chem.species.2.atoms", "1,8"); // H–O
    c.set("chem.species.2.bonds", "0-1:1");

    // One exothermic reaction. Enthalpy derived from bonds (negative = exothermic).
    // A low barrier and a modest pre-exponential so it fires at ~1000 K.
    c.set("chem.reaction.count", "1");
    c.set("chem.reaction.0.reactants", "0:1,1:1");
    c.set("chem.reaction.0.products", "2:2");
    c.set("chem.reaction.0.activation_j", "20000");
    c.set("chem.reaction.0.pre_exponential", "1.0e5");
    // enthalpy_j omitted → derived from the bond model.
    c
}

fn reacting_sim(seed: u64, temp_k: f64, n_reactant: i128) -> Simulation {
    let cfg = reacting_config(seed);
    let mut sim = Simulation::new(cfg);
    let spec = GridSpec::new(1, 1, 1, 1.0);
    let materials = sim.world.materials.clone();
    sim.world.set_generator(Box::new(GasParcelGen {
        spec,
        temp_k,
        n_reactant,
        n_species: 3,
        rho: 1000.0, // dense enough that thermal energy is well above the µJ quantum
        materials,
    }));
    sim.world.activate(ChunkCoord::new(0, 0, 0));
    sim
}

fn cell_energy(sim: &Simulation) -> i128 {
    sim.world
        .chunks
        .iter()
        .filter_map(|(_, c)| c.columns().get::<i128>(PHYS_ENERGY))
        .flat_map(|col| col.as_slice().iter().copied())
        .sum()
}

fn species_total(sim: &Simulation, s: usize) -> i128 {
    sim.world
        .chunks
        .iter()
        .filter_map(|(_, c)| c.columns().get::<i128>(species_column(s)))
        .flat_map(|col| col.as_slice().iter().copied())
        .sum()
}

/// **The Phase 3b headline.** An exothermic reaction, running inside the
/// simulation loop, raises the energy of the *shared* thermal field — the same
/// column physics reads. No external heat is added (flux off, radiation off), so
/// any energy gain came from bonds becoming heat.
#[test]
fn chemistry_heats_the_shared_thermal_field_in_the_sim_loop() {
    // Sized so the heat is *resolvable*, which is arithmetic rather than a fudge.
    // The thermal field is quantised to the microjoule and one bond is ~5e-19 J,
    // so fewer than about 2e12 reaction events cannot move it however long the
    // run. This test used 4e6 molecules and passed only because bond-derived
    // enthalpies were 6.022e23 times too large — it was measuring the units bug,
    // not the coupling. The coupling is real; it just needs a real number of
    // molecules, which is the honest statement about a microjoule-quantised
    // field talking to molecule-scale chemistry.
    let mut sim = reacting_sim(1, 1000.0, 40_000_000_000_000);

    let e0 = cell_energy(&sim);
    let product0 = species_total(&sim, 2);
    assert_eq!(product0, 0, "should start with no product");

    sim.run(200);

    let e1 = cell_energy(&sim);
    let product1 = species_total(&sim, 2);

    assert!(product1 > 0, "the reaction never fired inside the sim loop");
    assert!(
        e1 > e0,
        "chemistry released heat but the shared thermal field did not rise: {} → {} µJ",
        e0,
        e1
    );
}

/// Atoms are conserved through the coupled run. Hydrogen and oxygen are only
/// rearranged between H₂/O₂/H–O; none is created or destroyed by the reaction
/// running in the loop. (H atoms = 2·H₂ + 1·H–O; O atoms = 2·O₂ + 1·H–O.)
#[test]
fn atoms_are_conserved_through_the_coupled_run() {
    let mut sim = reacting_sim(2, 1000.0, 4_000_000);

    let h0 = 2 * species_total(&sim, 0) + species_total(&sim, 2);
    let o0 = 2 * species_total(&sim, 1) + species_total(&sim, 2);

    sim.run(300);

    let h1 = 2 * species_total(&sim, 0) + species_total(&sim, 2);
    let o1 = 2 * species_total(&sim, 1) + species_total(&sim, 2);

    assert_eq!(h0, h1, "hydrogen atoms were not conserved");
    assert_eq!(o0, o1, "oxygen atoms were not conserved");
}

/// Determinism survives chemistry. Two runs of the same seed and config produce
/// bit-identical worlds — the world hash, which now includes chemistry's events
/// and the species populations in every dirty chunk, is the same.
#[test]
fn chemistry_is_deterministic() {
    let run = || {
        let mut sim = reacting_sim(7, 1000.0, 4_000_000);
        sim.run(150);
        sim.world.world_hash()
    };
    assert_eq!(run(), run(), "the same seed produced two different chemical histories");
}

/// Save/resume is identical with chemistry active. A world saved mid-reaction,
/// reloaded, and finished must equal one that never stopped — including the
/// species populations and the shared energy the reaction has been changing. This
/// is the test that proves chemistry state round-trips through the snapshot codec.
#[test]
fn chemistry_survives_save_and_resume() {
    // Run to a midpoint and snapshot.
    let mut sim = reacting_sim(9, 1000.0, 4_000_000);
    sim.run(80);
    let bytes = sim.save();

    // Continue the original to the end.
    sim.run(80);
    let hash_continuous = sim.world.world_hash();

    // Resume from the snapshot and finish. The resumed world must reinstall the
    // generator (as the physics resume test does) so any chunk it needs to
    // regenerate matches — but our active chunk is dirty and lives in the
    // snapshot, so its chemistry comes back from disk, not from the generator.
    let mut resumed = Simulation::load(&bytes, reacting_config(9)).expect("resume must succeed");
    let spec = GridSpec::new(1, 1, 1, 1.0);
    let materials = resumed.world.materials.clone();
    resumed.world.set_generator(Box::new(GasParcelGen {
        spec,
        temp_k: 1000.0,
        n_reactant: 4_000_000,
        n_species: 3,
        rho: 1000.0,
        materials,
    }));
    resumed.run(80);
    let hash_resumed = resumed.world.world_hash();

    assert_eq!(
        hash_continuous, hash_resumed,
        "a world resumed mid-reaction diverged from one that never stopped"
    );
}

/// A world with chemistry declared but no reactants present does nothing and adds
/// nothing to history — a dormant chemistry is invisible, so worlds that never
/// react are unaffected by the layer existing. (This mirrors the physics
/// guarantee that an empty active set changes nothing.)
#[test]
fn a_chemistry_with_no_reactants_is_inert() {
    // Seed zero reactant molecules.
    let mut sim = reacting_sim(3, 1000.0, 0);
    let e0 = cell_energy(&sim);
    let hash0_events = sim.world.events.total();

    sim.run(100);

    assert_eq!(cell_energy(&sim), e0, "energy changed with nothing to react");
    assert_eq!(species_total(&sim, 2), 0, "product appeared from nothing");
    assert_eq!(
        sim.world.events.total(),
        hash0_events,
        "a dormant chemistry emitted events"
    );
}
