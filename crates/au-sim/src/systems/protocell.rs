//! The life-chemistry system — where bags nucleate, feed, divide and burst.
//!
//! Twin of `systems::chemistry` in shape: gather the bulk columns, compute,
//! scatter back. What is different is that the state it owns is not a field. It
//! is a list of individuals, and individuals are created and destroyed.
//!
//! # Order within a tick, and why
//!
//! ```text
//!   nucleate  ->  take up  ->  decide fate  ->  divide or burst
//! ```
//!
//! Uptake before fate, because a bag's fate is a statement about the bag *as it
//! now is*: deciding first and then feeding would let a vesicle divide on last
//! tick's geometry. Nucleation first, because a bag that closes this tick should
//! get to feed this tick — otherwise the moment of its birth is a tick in which
//! it does nothing, which is arbitrary.
//!
//! # Nothing here creates life
//!
//! There is a real risk in this file of writing `spawn_protocell()` and calling
//! it emergence. Every event below is a consequence of a quantity that already
//! existed:
//!
//! * **Nucleation** happens when a grid cell's own assembled surfactant — the
//!   quantity Phase 5f already computes, from the CMC, from amphiphilicity
//!   scores read off discovered graphs — exceeds what it takes to close the
//!   smallest bag a bilayer can bend into. The material comes out of the cell.
//!   Nothing is conjured; a surface that was already there closes.
//! * **Uptake** is Fick's law across a membrane of known area, thickness and
//!   permeability, driven by the concentration difference the bag itself
//!   creates.
//! * **Fate** is the geometry of step 1: `v > 1` bursts, `v <= 1/√2` divides.
//!   No size, no timer, no counter.
//! * **Division** partitions by a fair coin per molecule.
//!
//! The one declared quantity is the minimum vesicle radius, and it is declared
//! because the model has no bending rigidity to derive it from. It is set at
//! five bilayer thicknesses (~20 nm), which is about where real bilayers stop
//! being able to curve. Stated here rather than buried.

use au_core::schedule::{System, SystemDesc};
use au_core::{Domain, Layer, Rng};
use au_data::chunk::ChunkCoord;
use au_physics::GridSpec;

use au_chem::{react, CellChemistry, ReactingCell, Reaction, SpeciesId};

use crate::chem_columns::{has_chemistry, species_column};
use crate::protocell::Protocell;
use crate::world::World;

/// The dials of this layer. All but one are inherited from Phase 5f.
#[derive(Clone, Copy, Debug)]
pub struct ProtocellRules {
    pub membrane: au_chem::MembraneRules,
    /// Bilayer thickness, m. Real bilayers are about 4 nm across.
    pub thickness_m: f64,
    /// Areal strain a membrane survives before failing. ~3% for real lipids.
    pub lysis_strain: f64,
    /// Smallest bag a bilayer can bend into, as a multiple of its thickness.
    /// The one declared number in the file; see the module docs.
    pub min_radius_thicknesses: f64,
}

impl Default for ProtocellRules {
    fn default() -> Self {
        ProtocellRules {
            membrane: au_chem::MembraneRules::default(),
            thickness_m: 4.0e-9,
            lysis_strain: 0.03,
            min_radius_thicknesses: 5.0,
        }
    }
}

/// Volume one molecule occupies, m³ — roughly that of a water molecule. Sets the
/// densest a bag can physically be.
const MOLECULAR_VOLUME_M3: f64 = 3.0e-29;

/// Boltzmann's constant, J/K. Measured, not chosen.
const BOLTZMANN: f64 = 1.380_649e-23;

/// Bags one cell may close in one tick.
///
/// A **declared** quantity, and it bounds computation rather than physics: a
/// cell holding surface for a million vesicles would otherwise allocate a
/// million entities in one step. Surface left over closes on the next tick, so
/// the only thing this changes is how finely the work is spread.
const MAX_NUCLEATIONS_PER_CELL_TICK: i128 = 64;

pub struct ProtocellSystem {
    rules: ProtocellRules,
    n_species: usize,
    spec: GridSpec,
    cell_volume_m3: f64,
    /// D₀ of a 1-dalton reference particle — the same dial species transport uses.
    diffusion_d0: f64,
    /// Amphiphilicity per species, extended as the world invents molecules.
    species_amph: Vec<f64>,
    /// Molecular mass per species, for Graham scaling of the uptake rate.
    species_mass_mda: Vec<f64>,
    /// Scratch: bulk populations for one chunk.
    bulk: Vec<Vec<i128>>,
    /// The reaction set protocells run on their own contents.
    ///
    /// **This is the whole point of the step.** A bag that only absorbed and
    /// split would be a crystal with a membrane: its composition would be
    /// inherited and would *do* nothing, so variation could never affect
    /// persistence and no adaptation could occur. Chemistry inside the bag is
    /// what closes the loop — composition determines what the bag makes,
    /// what it makes determines how fast it divides, and division propagates
    /// the composition.
    ///
    /// Declared reactions only, for now. Open-ended generation *inside*
    /// protocells needs the reaction cache shared between two systems rather
    /// than rebuilt in each, and that is a real refactor rather than a line.
    /// Stated as a limit, not hidden.
    reactions: Vec<Reaction>,
    /// Scratch: one protocell's contents as a dense population vector.
    inner: CellChemistry,
    /// Cells that are open ports to an infinitely dilute exterior (Phase 5h).
    ///
    /// **A drain removes whatever is in it, and a bag is a thing in it.** Until
    /// now protocells could only die by bursting, so a population had births and
    /// no removal — it could only accumulate. That is fatal to selection: with no
    /// turnover, a lineage that divides faster never *displaces* anything, it
    /// merely adds to a pile. Washout is what makes reproduction rate matter.
    ///
    /// The drain stays blind, exactly as it is for molecules. It removes bags
    /// without regard to what they are made of, which is what keeps selection an
    /// outcome rather than an input.
    sink_cells: Vec<usize>,
    /// Atoms per species, for booking what leaves.
    species_formula: Vec<Vec<i128>>,
    /// Dynamic viscosity of the medium, Pa·s. Read from the material table at
    /// boot; a protocell's Brownian mobility is inversely proportional to it.
    medium_viscosity_pa_s: f64,
    /// Scratch: cell velocity, m/s per axis.
    vel: Vec<[f64; 3]>,

    pub last_drifted: u64,
    pub last_reacted: u64,
    pub last_nucleated: u64,
    pub last_divided: u64,
    pub last_lysed: u64,
}

impl ProtocellSystem {
    pub fn from_config(world: &World) -> Option<(SystemDesc, ProtocellSystem)> {
        if !world.config.bool_or("life.protocells", false) {
            return None;
        }
        let chem = crate::chemistry::Chemistry::from_kv(world.config.as_map()).ok()??;
        let d = ProtocellRules::default();
        let rules = ProtocellRules {
            membrane: chem.membrane.unwrap_or_default(),
            thickness_m: world.config.f64_or("life.bilayer_thickness_m", d.thickness_m),
            lysis_strain: world.config.f64_or("life.lysis_strain", d.lysis_strain),
            min_radius_thicknesses: world
                .config
                .f64_or("life.min_radius_thicknesses", d.min_radius_thicknesses),
        };
        let spec = GridSpec::new(
            world.config.u64_or("physics.grid.nx", 1) as usize,
            world.config.u64_or("physics.grid.ny", 1) as usize,
            world.config.u64_or("physics.grid.nz", 1) as usize,
            world.config.f64_or("physics.grid.cell_m", 1.0),
        );
        let n_species = chem.column_count();

        // Life chemistry sits above Chemistry in the stack, so it may read
        // downward into the species fields — and must run after them, which the
        // layer declaration is what guarantees.
        let desc = SystemDesc::new("life.protocells", Layer::LifeChemistry)
            .reads(&[Layer::Universe, Layer::Physics, Layer::Chemistry, Layer::LifeChemistry])
            .every(1);

        Some((
            desc,
            ProtocellSystem {
                rules,
                n_species,
                spec,
                cell_volume_m3: chem.cell_volume_m3,
                diffusion_d0: chem.diffusion_m2_s,
                species_amph: Vec::new(),
                species_mass_mda: Vec::new(),
                vel: Vec::new(),
                // Pa·s already — `mu_liquid` is dynamic viscosity, as its own
                // documentation says. Read from the same material table
                // everything else uses; a world whose medium has no declared
                // viscosity gets no Brownian motion, which is the correct
                // reading of "nobody said what this fluid is like".
                medium_viscosity_pa_s: world
                    .materials
                    .by_name("water")
                    .and_then(|id| world.materials.get(id))
                    .and_then(|m| m.viscosity(au_physics::Phase::Liquid))
                    .unwrap_or(0.0),
                reactions: chem.reactions.clone(),
                inner: CellChemistry::new(n_species),
                sink_cells: chem
                    .ports
                    .iter()
                    .filter(|(_, p)| matches!(p, au_chem::Port::Sink))
                    .map(|(c, _)| *c)
                    .collect(),
                species_formula: Vec::new(),
                bulk: vec![Vec::new(); n_species],
                last_drifted: 0,
                last_reacted: 0,
                last_nucleated: 0,
                last_divided: 0,
                last_lysed: 0,
            },
        ))
    }

    /// Molecules needed to close a sphere of radius `r`.
    fn molecules_for_radius(&self, r: f64) -> i128 {
        let area = 4.0 * std::f64::consts::PI * r * r;
        (area / self.rules.membrane.area_per_molecule_m2).ceil() as i128
    }

    /// Is this species part of a membrane rather than a solute?
    fn is_membrane(&self, s: usize) -> bool {
        self.species_amph.get(s).copied().unwrap_or(0.0) >= self.rules.membrane.min_score
    }

    /// A protocell's membrane count, osmolyte count, and shape.
    fn shape_of(&self, p: &Protocell, c_external: f64) -> Option<au_chem::Shape> {
        let mut membrane = 0i128;
        let mut osmolytes = 0i128;
        for (s, n) in &p.contents {
            if self.is_membrane(*s as usize) {
                membrane += *n;
            } else {
                osmolytes += *n;
            }
        }
        au_chem::shape(membrane, osmolytes, c_external, self.rules.membrane.area_per_molecule_m2)
    }
}

impl System<World> for ProtocellSystem {
    fn run(&mut self, world: &mut World) {
        if self.n_species == 0 {
            return;
        }
        let dt = world.clock.scale().as_secs_f64();
        if dt <= 0.0 {
            return;
        }
        let tick = world.clock.tick().0;
        let seed = world.seed;

        // Per-species derived properties, extended as the vocabulary grows. A
        // molecule cannot change once interned, so old entries stay right.
        let known = world.chem_registry.len().min(self.n_species);
        while self.species_amph.len() < known {
            let id = au_chem::SpeciesId(self.species_amph.len() as u32);
            let m = world.chem_registry.get(id).unwrap();
            let table = {
                // The table lives in the chemistry description; rebuilding it per
                // species would be absurd, so it is fetched once per growth step.
                match crate::chemistry::Chemistry::from_kv(world.config.as_map()) {
                    Ok(Some(c)) => c.table,
                    _ => break,
                }
            };
            self.species_amph.push(au_chem::amphiphilicity(m, &table));
            self.species_mass_mda.push(m.mass_mda(&table).max(1) as f64);
        }

        if std::env::var("AU_DBG").is_ok() && tick == 1 {
            for (s, sc) in self.species_amph.iter().enumerate() {
                eprintln!(
                    "SPECIES {} amph={:.4} min={:.4} membrane={}",
                    s, sc, self.rules.membrane.min_score, *sc >= self.rules.membrane.min_score
                );
            }
        }
        if !self.sink_cells.is_empty() {
            let table = match crate::chemistry::Chemistry::from_kv(world.config.as_map()) {
                Ok(Some(c)) => Some(c.table),
                _ => None,
            };
            if let Some(t) = table {
                if self.species_formula.len() < world.chem_registry.len() {
                    self.species_formula = au_chem::formula_table(&world.chem_registry, &t);
                }
            }
        }

        let coords: Vec<ChunkCoord> = world.active.iter().copied().collect();
        let mut nucleated = 0u64;
        let mut divided = 0u64;
        let mut lysed = 0u64;
        let mut washed = 0u64;
        let mut drifted = 0u64;
        let mut reacted = 0u64;
        let mut heat_out: i128 = 0;

        // Protocells are taken out of the store, worked on, and put back. Doing
        // it in place would mean holding a borrow of the store while allocating
        // ids from the world, and the borrow checker is right to refuse: a
        // division mutates both at once.
        let mut store = std::mem::take(&mut world.protocells);

        for coord in &coords {
            let Some(chunk) = world.chunks.get_mut(*coord) else { continue };
            if !has_chemistry(chunk.columns(), self.n_species) {
                continue;
            }
            let cells = self.spec.cells();

            // ── Gather the bulk. ────────────────────────────────────────────
            {
                let cols = chunk.columns();
                for s in 0..self.n_species {
                    self.bulk[s].clear();
                    self.bulk[s]
                        .extend_from_slice(cols.get::<i128>(species_column(s)).unwrap().as_slice());
                }
            }

            // Temperature per cell, derived exactly as physics derives it, so the
            // two layers never disagree about how hot anything is.
            let mut temp_k = vec![0.0f64; cells];
            {
                let cols = chunk.columns();
                let mass = cols.get::<i128>(crate::physics_columns::PHYS_MASS).unwrap().as_slice();
                let energy =
                    cols.get::<i128>(crate::physics_columns::PHYS_ENERGY).unwrap().as_slice();
                let matid =
                    cols.get::<u16>(crate::physics_columns::PHYS_MATERIAL).unwrap().as_slice();
                for c in 0..cells {
                    let kg = mass[c] as f64 / au_physics::AG_PER_KG as f64;
                    if kg <= 0.0 {
                        continue;
                    }
                    if let Some(m) = world.materials.get(au_physics::MaterialId(matid[c])) {
                        let st = au_physics::derive(
                            kg,
                            au_physics::Energy(energy[c]).as_joules(),
                            m,
                            0.0,
                            self.cell_volume_m3,
                        );
                        temp_k[c] = st.temperature;
                        if std::env::var("AU_DBG").is_ok() && tick == 1 {
                            eprintln!(
                                "TEMP raw_mass={} kg={:.4e} raw_e={} e_j={:.4e} matid={} T={:.4e} cap={:.4e}",
                                mass[c], kg, energy[c],
                                au_physics::Energy(energy[c]).as_joules(),
                                matid[c], st.temperature, st.heat_capacity
                            );
                        }
                    }
                }
            }

            // ── The flow. ───────────────────────────────────────────────────
            //
            // Velocity is derived from exact momentum and exact mass, never
            // stored — the same discipline as temperature. A chunk with no fluid
            // solver has no momentum columns and every bag stays where it is,
            // which is the correct reading of a world with no currents in it.
            self.vel.clear();
            self.vel.resize(cells, [0.0; 3]);
            {
                let cols = chunk.columns();
                let mass = cols.get::<i128>(crate::physics_columns::PHYS_MASS).unwrap().as_slice();
                for (k, id) in [
                    crate::physics_columns::PHYS_MOM_X,
                    crate::physics_columns::PHYS_MOM_Y,
                    crate::physics_columns::PHYS_MOM_Z,
                ]
                .iter()
                .enumerate()
                {
                    if let Some(col) = cols.get::<i128>(*id) {
                        let m = col.as_slice();
                        for c in 0..cells.min(m.len()).min(mass.len()) {
                            let kg = mass[c] as f64 / au_physics::AG_PER_KG as f64;
                            if kg > 0.0 {
                                self.vel[c][k] = au_physics::Momentum(m[c]).as_si() / kg;
                            }
                        }
                    }
                }
            }

            // External concentration of solute, per cell. Membrane species are
            // excluded: a surfactant sitting in a surface is not osmotically
            // active, which is the same reason Phase 5g excluded it from
            // diffusion.
            let mut c_ext = vec![0.0f64; cells];
            for (cell, c) in c_ext.iter_mut().enumerate() {
                let mut solute = 0i128;
                for s in 0..self.n_species {
                    if !self.is_membrane(s) {
                        solute += self.bulk[s][cell];
                    }
                }
                *c = (solute.max(0) as f64) / self.cell_volume_m3;
            }

            // ── Nucleation. ─────────────────────────────────────────────────
            //
            // A cell whose assembled surfactant exceeds what it takes to close
            // the smallest bendable bag closes one. The material is taken out of
            // the cell; the bag captures ambient solute in proportion to the
            // volume it encloses. Exactly conservative in both directions.
            let r_min = self.rules.min_radius_thicknesses * self.rules.thickness_m;
            let n_min = self.molecules_for_radius(r_min);
            for cell in 0..cells {
                let mut assembled = 0i128;
                for s in 0..self.n_species {
                    if !self.is_membrane(s) {
                        continue;
                    }
                    assembled += au_chem::assembled_count(
                        self.bulk[s][cell],
                        self.species_amph[s.min(self.species_amph.len().saturating_sub(1))],
                        self.cell_volume_m3,
                        &self.rules.membrane,
                    );
                }
                if assembled < n_min {
                    continue;
                }
                // ── How many bags close, not how many may. ──────────────
                //
                // The assembled surfactant in a cell is *already structure* —
                // Phase 5f caps free monomer at the CMC and everything above it
                // is surface, so 42,000,000 assembled molecules are 42,000,000
                // molecules of micelle and sheet, not a supersaturated solution
                // waiting to precipitate. (That was the first diagnosis, and
                // reading `assembled_count` disproved it: the model was right.)
                //
                // What was wrong was here. This loop closed *one* bag per cell
                // per tick, which is not a physical rate but the shape of a for
                // loop — and it left 3,300 bags' worth of surface lying around
                // as unclosed aggregate for thousands of ticks, so the medium
                // held the material its own compartments were made of.
                //
                // A sheet of amphiphile closes because a closed vesicle has no
                // exposed edge, and the edge is what costs energy. Nothing about
                // that says one per second. So every bag the surface can support
                // closes, bounded only by the material present.
                let mut budget = assembled;
                // ── How big is one vesicle? ─────────────────────────────
                //
                // The first version closed *all* the assembled surfactant into a
                // single bag, and the result was a vesicle eleven times larger
                // than the cell containing it: it swallowed the entire medium,
                // fissioned, and its daughters then burst because the emptied
                // medium made their osmotic volume explode. Everything dumped
                // back, it renucleated, and the world sat in a two-tick cycle
                // that no chemistry could influence.
                //
                // The error was physical rather than numerical. Surfactant above
                // the CMC does not form one enormous vesicle; it forms **many
                // small ones**, because the bending energy of a nearly-flat
                // membrane is what would have to pay for a large one. So one
                // nucleation takes exactly the molecules needed to close the
                // smallest bag a bilayer can bend into, and leaves the rest in
                // the medium to close bags of their own on later ticks.
                //
                // A bag born this way is a sphere — `v = 1` — and therefore
                // stable. It can only reach the fission bound by *growing*, and
                // since membrane material cannot cross a membrane, the only way
                // it can grow is by making surfactant itself. That is the loop.
                let area = n_min as f64 * self.rules.membrane.area_per_molecule_m2;
                let r = (area / (4.0 * std::f64::consts::PI)).sqrt();
                let v_bag = (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
                let frac = (v_bag / self.cell_volume_m3).min(1.0);

                // Cap per tick so one cell cannot allocate a million entities in
                // a single step. Declared, and stated as declared: it bounds the
                // *work*, not the physics, and a cell with more surface than this
                // simply keeps closing bags on later ticks.
                let per_tick = budget / n_min.max(1);
                let per_tick = per_tick.min(MAX_NUCLEATIONS_PER_CELL_TICK);

                for _ in 0..per_tick {
                if budget < n_min {
                    break;
                }
                budget -= n_min;
                let mut p = Protocell {
                    id: world.entities.alloc(),
                    parent: au_core::ids::EntityId::NONE,
                    born_tick: tick,
                    chunk: *coord,
                    cell: cell as u32,
                    contents: Vec::new(),
                };
                // ⚠ MEASURED DEFECT, NOT YET FIXED. Five facts. Read all of them
                // before touching anything here.
                //
                // 1. **Nucleation itself is correct.** Traced: captured 1675
                //    osmolytes against an osmotic balance of 1675.6 — ratio
                //    1.000. A bag is born a sphere at `v ≈ 1`, as intended.
                //
                // 2. **Fissioning bags have wandered.** They nucleate in the cell
                //    where surfactant is seeded and fission seven or eight cells
                //    away. Brownian motion is implicated.
                //
                // 3. **They lose a third of their solute in transit.** `osmo`
                //    1675 at birth, 1132 at fission, `memb` still exactly 12567
                //    and the area unchanged. Losing solute at fixed area is what
                //    drops `v` below 1/√2.
                //
                // 4. **No chemistry is involved.** 46,497 divisions from 1,020
                //    nucleations in a world with *no reactions declared*.
                //
                // 5. **The mechanism.** Uptake removes it, and `c_in == c_ext`
                //    fails for a collective reason. Traced:
                //
                //        FLUX sp=2 c_in=5.00000e25 c_out=4.98929e25 n=-3
                //
                //    Each nucleation takes ~1675 osmolytes out of the medium, so
                //    `c_ext` falls; every bag already floating there is then above
                //    the new external concentration and sheds solute to match.
                //    Its osmotic volume shrinks at fixed area, `v` crosses 1/√2,
                //    and it fissions. **Nucleation triggers a fission wave across
                //    the standing population**, and each wave doubles it.
                //
                //    Part is real physics — draw solute from a medium and vesicles
                //    shrink. What is not real is the *coupling*: bags are made from
                //    the same finite cell they float in, so their own creation
                //    squeezes their neighbours. A real medium is not measurably
                //    depleted by one vesicle forming in it. Candidate fixes: a
                //    reservoir holding solute fixed (as Phase 5h holds substrate),
                //    or nucleation drawing osmolytes from outside the cell budget.
                //    Not chosen, because choosing without measuring produced four
                //    wrong hypotheses before this one.
                //
                // **Until this is closed, division counts are not evidence of
                // growth in any world**, and the Phase 5i step 2 claim that
                // composition determines replication rate must be re-earned
                // against it rather than assumed to have survived.
                //
                // Surfactant: `n_min` molecules, drawn from the assembled pool in
                // species order. Deterministic, and it prefers whichever species
                // the world declared first — arbitrary, but arbitrary in a way
                // that is stated rather than random.
                let mut need = n_min;
                for s in 0..self.n_species {
                    let take = if self.is_membrane(s) {
                        let avail = au_chem::assembled_count(
                            self.bulk[s][cell],
                            self.species_amph[s.min(self.species_amph.len().saturating_sub(1))],
                            self.cell_volume_m3,
                            &self.rules.membrane,
                        );
                        let t = avail.min(need);
                        need -= t;
                        t
                    } else {
                        // Solute: whatever was in the volume the bag closed around.
                        ((self.bulk[s][cell] as f64) * frac) as i128
                    };
                    let take = take.min(self.bulk[s][cell]).max(0);
                    if take > 0 {
                        self.bulk[s][cell] -= take;
                        p.add(s as u32, take);
                    }
                }
                if std::env::var("AU_DBG").is_ok() {
                    eprintln!(
                        "NUC r={:.3e} v_bag={:.3e} frac={:.3e} took={}",
                        r, v_bag, frac, p.total()
                    );
                }
                store.insert(p);
                nucleated += 1;
                }
            }

            // ── Uptake, then fate. ──────────────────────────────────────────
            let ids: Vec<_> = store
                .iter()
                .filter(|p| p.chunk == *coord)
                .map(|p| p.id)
                .collect();
            for id in ids {
                let Some(p) = store.get(id).cloned() else { continue };
                let cell = p.cell as usize;
                if cell >= cells {
                    continue;
                }
                let Some(sh) = self.shape_of(&p, c_ext[cell].max(1.0)) else { continue };

                // Permeability of this bag's own surface. A bag is by definition
                // closed, so its coverage is complete and its permeability is the
                // residual leak of a real bilayer — which is what lets it take
                // anything up at all.
                let perm = au_chem::permeability(1.0, &self.rules.membrane);

                // Fick across the membrane: J = (D·p/d)·A·Δc, in molecules.
                let mut updated = p.clone();
                if self.diffusion_d0 > 0.0 && sh.volume_m3 > 0.0 {
                    for s in 0..self.n_species {
                        // Membrane species are structure, not solute. Phase 5g
                        // made exactly this exclusion for grid transport — a
                        // molecule that is part of a surface has no independent
                        // gradient to run down — and the same holds here. Without
                        // it a bag's own surfactant leaks into a medium sitting at
                        // the CMC and the vesicle dissolves itself within a few
                        // ticks of forming.
                        //
                        // The consequence is deliberate: a bag cannot grow its
                        // membrane by absorbing monomer, so **the only way it can
                        // grow is by making surfactant itself**. That is the loop
                        // this step exists to close, and it is better that it be
                        // the only route than one of two.
                        if self.is_membrane(s) {
                            continue;
                        }
                        let mass = self.species_mass_mda.get(s).copied().unwrap_or(1000.0);
                        let d_s = self.diffusion_d0 * (1000.0 / mass).sqrt();
                        let c_in = updated.count(s as u32) as f64 / sh.volume_m3;
                        let c_out = self.bulk[s][cell] as f64 / self.cell_volume_m3;
                        let flux =
                            (d_s * perm / self.rules.thickness_m) * sh.area_m2 * (c_out - c_in) * dt;
                        // ── No face may overshoot equilibrium. ──────────────
                        //
                        // This is an explicit Euler step, and an explicit step
                        // across a membrane is stable only while `dt` is shorter
                        // than the bag's equilibration time — here about 0.15 s
                        // against a tick of 1 s. Without a limiter the flux
                        // sails past equality and reverses next tick, and a
                        // daughter whose parent's chemistry ate all its substrate
                        // sees a near-infinite gradient and floods: 1,563
                        // molecules into an osmolyte budget of 837, which bursts
                        // it. Every bag in the world was dying of this.
                        //
                        // The cure is the one this engine has used three times
                        // already — Phase 2's CFL substepping, Phase 5b's
                        // donor's-pocket limiter, Phase 4b's convex combination:
                        // **transport may be wrong, but it may not be
                        // impossible.** Cap the transfer at the amount that
                        // brings the two sides to equal concentration. Two
                        // well-mixed volumes in contact equalise after
                        // `Δc · V_in · V_out / (V_in + V_out)` molecules, and the
                        // bag is vanishingly small against its cell, so that is
                        // `Δc · V_in` to a part in a million.
                        let equalising =
                            ((c_out - c_in).abs() * sh.volume_m3 * self.cell_volume_m3
                                / (sh.volume_m3 + self.cell_volume_m3))
                                .floor();
                        let mut n = flux.abs().min(equalising).copysign(flux) as i128;
                        if n > 0 {
                            n = n.min(self.bulk[s][cell]);
                        } else if n < 0 {
                            n = n.max(-updated.count(s as u32));
                        }
                        if n != 0 {
                            self.bulk[s][cell] -= n;
                            updated.add(s as u32, n);
                        }
                    }
                }

                // ── Chemistry inside the bag. ───────────────────────────────
                //
                // At the bag's *own* volume, which is the physically important
                // part: a vesicle is small, so its contents are concentrated, so
                // bimolecular rates inside it are far higher than in the medium
                // it floats in. Concentration is the classical reason
                // compartments matter for the origin of life, and here it is not
                // asserted — it falls out of dividing by a smaller volume.
                //
                // Temperature comes from the medium: a bilayer is four
                // nanometres of hydrocarbon and does not insulate. Reaction heat
                // goes back into the grid cell's thermal energy for the same
                // reason, so energy stays exactly conserved across the boundary.
                // A bag cannot be denser than its own molecules allow. Without
                // this floor a nearly-empty vesicle reports a near-zero osmotic
                // volume, its contents look infinitely concentrated, and the
                // reactor — correctly — subdivides the timestep towards infinity
                // trying to hold every reaction under its consumption cap. The
                // first run of these tests hung on exactly that.
                //
                // Excluded volume is real physics, not a numerical patch:
                // molecules occupy space, and a bag of N of them is at least
                // N·v_molecule across.
                let packed = updated.total().max(0) as f64 * MOLECULAR_VOLUME_M3;
                let react_volume = sh.volume_m3.max(packed);
                if std::env::var("AU_DBG").is_ok() && tick % 200 == 1 {
                    eprintln!(
                        "GATE n_rxn={} vol={:.3e} T={:.2}",
                        self.reactions.len(), react_volume, temp_k[cell]
                    );
                }
                if !self.reactions.is_empty() && react_volume > 0.0 && temp_k[cell] > 0.0 {
                    self.inner.ensure_species(self.n_species);
                    for s in 0..self.n_species {
                        self.inner.set(SpeciesId(s as u32), updated.count(s as u32));
                    }
                    let mut heat = au_physics::Energy(0);
                    let mut rc = ReactingCell {
                        chem: &mut self.inner,
                        thermal_energy: &mut heat,
                        temperature_k: temp_k[cell],
                        volume_m3: react_volume,
                    };
                    let report = react(&mut rc, &self.reactions, dt);
                    if std::env::var("AU_DBG").is_ok() && tick % 200 == 1 {
                        eprintln!(
                            "RXN n_rxn={} T={:.1} V={:.3e} dt={} events={} inner0={} inner1={} inner2={}",
                            self.reactions.len(), temp_k[cell], react_volume, dt, report.events,
                            self.inner.get(SpeciesId(0)), self.inner.get(SpeciesId(1)),
                            self.inner.get(SpeciesId(2))
                        );
                    }
                    if report.events > 0 {
                        updated.contents.clear();
                        for s in 0..self.n_species {
                            let n = self.inner.get(SpeciesId(s as u32));
                            if n != 0 {
                                updated.add(s as u32, n);
                            }
                        }
                        reacted += report.events;
                        heat_out += heat.0;
                    }
                }

                let Some(sh) = self.shape_of(&updated, c_ext[cell].max(1.0)) else {
                    store.insert(updated);
                    continue;
                };
                if std::env::var("AU_DBG").is_ok() {
                    eprintln!(
                        "BAG id={} v={:.4} area={:.3e} vol={:.3e} c_ext={:.3e} total={}",
                        id.index, sh.reduced_volume, sh.area_m2, sh.volume_m3,
                        c_ext[cell], updated.total()
                    );
                }
                match au_chem::fate(sh.reduced_volume, self.rules.lysis_strain) {
                    au_chem::Fate::Intact => {
                        store.insert(updated);
                    }
                    au_chem::Fate::Lysis => {
                        // The bag fails and gives everything back to the medium.
                        for (s, n) in &updated.contents {
                            self.bulk[*s as usize][cell] += *n;
                        }
                        store.remove(id);
                        world.entities.free(id);
                        lysed += 1;
                    }
                    au_chem::Fate::Fission => {
                        // Enough surface for two spheres. A fair coin per
                        // molecule, from a stream derived from this bag's own id
                        // and this tick — so a resumed world splits identically.
                        let mut rng = Rng::derive(
                            seed,
                            Domain::CHEMISTRY,
                            id.index as u64 | ((id.generation as u64) << 32),
                            tick,
                        );
                        let counts: Vec<i128> = updated.contents.iter().map(|(_, n)| *n).collect();
                        let (a, b) = au_chem::partition(&counts, 0.5, &mut rng);

                        store.remove(id);
                        world.entities.free(id);
                        for half in [a, b] {
                            let mut child = Protocell {
                                id: world.entities.alloc(),
                                parent: id,
                                born_tick: tick,
                                chunk: *coord,
                                cell: cell as u32,
                                contents: Vec::new(),
                            };
                            for (i, (s, _)) in updated.contents.iter().enumerate() {
                                child.add(*s, half[i]);
                            }
                            store.insert(child);
                        }
                        divided += 1;
                    }
                }
            }

            // ── Bags wander. ────────────────────────────────────────────────
            //
            // **The transport that actually matters at this scale**, and it took
            // three sessions of chasing velocity fields to notice. A vesicle is a
            // colloid, and a colloid in a fluid is kicked by the molecules around
            // it. Stokes–Einstein gives the rate:
            //
            //     D = k_B·T / (6·π·η·r)
            //
            // For a 20 nm bag in water at 320 K that is 2.1e-11 m²/s, and the rms
            // displacement in one second is 6.5 µm — **six and a half cells of a
            // micron grid, per tick.** Advection by a convection roll is a rounding
            // error beside it. This is exactly why a bacterium needs a flagellum
            // to go anywhere in particular: below a few microns, swimming loses to
            // being shoved.
            //
            // Nothing here is declared. Temperature comes from the same `derive`
            // physics uses, viscosity from the material table, and the radius from
            // the bag's own membrane area — `r = √(A/4π)`, the area it has because
            // of how much surfactant it made. A bag that grows slows down, and it
            // slows as `1/r`, because that is what Stokes said.
            //
            // Same lattice treatment as advection below: a hop probability, so the
            // mean displacement is right and no fractional cell offset has to be
            // carried as accumulating state. Same derived RNG, so a resumed world
            // wanders identically.
            if dt > 0.0 && self.spec.cell_m > 0.0 {
                let wandering: Vec<_> = store
                    .iter()
                    .filter(|p| p.chunk == *coord)
                    .map(|p| (p.id, p.cell as usize))
                    .collect();
                for (id, cell) in wandering {
                    if cell >= cells || temp_k[cell] <= 0.0 {
                        continue;
                    }
                    let Some(p) = store.get(id).cloned() else { continue };
                    let Some(sh) = self.shape_of(&p, c_ext[cell].max(1.0)) else { continue };
                    let r_bag = (sh.area_m2 / (4.0 * std::f64::consts::PI)).sqrt();
                    if r_bag <= 0.0 {
                        continue;
                    }
                    // Dynamic viscosity from the host material's kinematic value.
                    let eta = self.medium_viscosity_pa_s;
                    if eta <= 0.0 {
                        continue;
                    }
                    let d_bag = BOLTZMANN * temp_k[cell]
                        / (6.0 * std::f64::consts::PI * eta * r_bag);
                    // One-dimensional rms step per axis, as a fraction of a cell.
                    let sigma = (2.0 * d_bag * dt).sqrt() / self.spec.cell_m;
                    if sigma <= 0.0 {
                        continue;
                    }
                    let mut rng = Rng::derive(
                        seed,
                        Domain::CHEMISTRY,
                        (id.index as u64) | ((id.generation as u64) << 32),
                        tick ^ 0xb2ff_00d5,
                    );
                    let (nx, ny) = (self.spec.nx, self.spec.ny);
                    let (mut cx, mut cy, mut cz) =
                        (cell % nx, (cell / nx) % ny, cell / (nx * ny));
                    let dims = [self.spec.nx, self.spec.ny, self.spec.nz];
                    for (axis, coordinate) in
                        [&mut cx, &mut cy, &mut cz].into_iter().enumerate()
                    {
                        if dims[axis] < 2 {
                            continue;
                        }
                        // A step of ±1 cell with probability σ²/2 each way keeps
                        // the variance right; beyond σ ≈ 1 the bag is crossing
                        // more than a cell per tick and the walk saturates, which
                        // is a resolution limit and is stated rather than hidden.
                        let p_side = (sigma * sigma / 2.0).min(0.5);
                        if rng.chance(p_side) {
                            if *coordinate + 1 < dims[axis] {
                                *coordinate += 1;
                            }
                        } else if rng.chance(p_side / (1.0 - p_side).max(1e-12)) {
                            if *coordinate > 0 {
                                *coordinate -= 1;
                            }
                        }
                    }
                    let dest = self.spec.idx(cx, cy, cz);
                    if dest != cell && dest < cells {
                        let mut moved = p;
                        moved.cell = dest as u32;
                        store.insert(moved);
                        drifted += 1;
                    }
                }
            }

            // ── Bags ride the fluid. ────────────────────────────────────────
            //
            // A vesicle is neutrally buoyant, so it goes where the water goes.
            // Over one tick it is displaced by `v·dt`, which is normally a small
            // fraction of a cell — and a lattice cannot represent a fractional
            // move.
            //
            // The honest treatment of a discrete particle on a lattice is a *hop
            // probability*: `p = |v|·dt/dx` toward the neighbour the flow points
            // at. The mean displacement is then exactly `v·dt`, which is the
            // analytic property worth having, and the variance is real dispersion
            // rather than an artefact — a cloud of bags in a shear does spread.
            //
            // The alternative was a sub-cell offset carried per bag, and it was
            // rejected: that is an accumulating `f64` in world state, and this
            // engine has kept floats to derived quantities since Phase 2 for
            // exactly the reason that accumulated float error is divergence. A
            // cell index stays an integer.
            //
            // Stochastic, deterministic, and not in tension: the draw comes from
            // `f(seed, domain, bag, tick)` as every other draw in the engine
            // does, so a resumed world moves its bags identically.
            if dt > 0.0 && self.spec.cell_m > 0.0 {
                let moving: Vec<_> = store
                    .iter()
                    .filter(|p| p.chunk == *coord)
                    .map(|p| (p.id, p.cell as usize))
                    .collect();
                for (id, cell) in moving {
                    if cell >= cells {
                        continue;
                    }
                    let v = self.vel[cell];
                    if v == [0.0; 3] {
                        continue;
                    }
                    let (nx, ny) = (self.spec.nx, self.spec.ny);
                    let (mut cx, mut cy, mut cz) = (cell % nx, (cell / nx) % ny, cell / (nx * ny));
                    let mut rng = Rng::derive(
                        seed,
                        Domain::CHEMISTRY,
                        (id.index as u64) | ((id.generation as u64) << 32),
                        tick ^ 0x5f5f_5f5f,
                    );
                    for (axis, coordinate) in [&mut cx, &mut cy, &mut cz].into_iter().enumerate() {
                        let f = v[axis] * dt / self.spec.cell_m;
                        if f == 0.0 {
                            continue;
                        }
                        if rng.chance(f.abs().min(1.0)) {
                            let n = [self.spec.nx, self.spec.ny, self.spec.nz][axis];
                            if f > 0.0 && *coordinate + 1 < n {
                                *coordinate += 1;
                            } else if f < 0.0 && *coordinate > 0 {
                                *coordinate -= 1;
                            }
                        }
                    }
                    let dest = self.spec.idx(cx, cy, cz);
                    if dest != cell && dest < cells {
                        if let Some(p) = store.get(id).cloned() {
                            let mut moved = p;
                            moved.cell = dest as u32;
                            store.insert(moved);
                            drifted += 1;
                        }
                    }
                }
            }

            // ── Washout. ────────────────────────────────────────────────────
            //
            // A bag sitting in an open port is gone, with everything it holds.
            // Its contents leave the world rather than returning to the medium —
            // that is the difference between washing out and bursting, and the
            // atom ledger records it so the books still balance exactly.
            if !self.sink_cells.is_empty() {
                let doomed: Vec<_> = store
                    .iter()
                    .filter(|p| p.chunk == *coord && self.sink_cells.contains(&(p.cell as usize)))
                    .map(|p| p.id)
                    .collect();
                for id in doomed {
                    if let Some(p) = store.remove(id) {
                        let mut lost = vec![0i128; self.n_species];
                        for (sp, n) in &p.contents {
                            if let Some(slot) = lost.get_mut(*sp as usize) {
                                *slot += *n;
                            }
                        }
                        world.atoms.book_gross(&[], &lost, &self.species_formula);
                        world.entities.free(id);
                        washed += 1;
                    }
                }
            }

            // ── Scatter the bulk back, and hand back the reaction heat. ────
            {
                let cols = chunk.columns_mut();
                if heat_out != 0 {
                    let e = cols.get_mut::<i128>(crate::physics_columns::PHYS_ENERGY).unwrap();
                    // Spread over the chunk's first cell that has mass: the bags
                    // are all in this chunk, and a bilayer conducts. Crude, and
                    // exact — no joule is created or lost.
                    e.as_mut_slice()[0] += heat_out;
                    heat_out = 0;
                }
                for s in 0..self.n_species {
                    cols.get_mut::<i128>(species_column(s))
                        .unwrap()
                        .as_mut_slice()
                        .copy_from_slice(&self.bulk[s]);
                }
            }
        }

        store.stat_nucleated += nucleated;
        store.stat_divided += divided;
        store.stat_lysed += lysed;
        store.stat_washed += washed;
        world.protocells = store;

        self.last_drifted = drifted;
        self.last_reacted = reacted;
        self.last_nucleated = nucleated;
        self.last_divided = divided;
        self.last_lysed = lysed;
    }
}
