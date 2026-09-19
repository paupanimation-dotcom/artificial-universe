//! The chemistry system — where reactions run in the real simulation loop.
//!
//! This is the twin of `systems::physics`, and it is deliberately built to the
//! same shape: gather columns into scratch, compute, scatter back, fold events and
//! the ledger into the world. Reading that file first makes this one obvious.
//!
//! # The one thing that makes Phase 3 real
//!
//! Chemistry does not have its own energy. It reads and writes the *same*
//! `PHYS_ENERGY` column that the thermodynamics solver conducts and radiates. So
//! when an exothermic reaction fires in a cell, the heat it releases lands in the
//! field physics will carry away on the next tick — and the temperature that
//! results feeds back into the Arrhenius rates on the tick after that. Heat drives
//! chemistry; chemistry drives heat; through one shared, conserved account.
//!
//! To make that coupling exact, the temperature chemistry uses for its rates is
//! derived from the shared energy by the *same* function physics uses
//! (`au_physics::derive`). Chemistry and physics therefore never disagree about
//! what temperature a given energy means — there is only one answer, computed one
//! way.
//!
//! # Ordering, and why it does not threaten determinism
//!
//! The scheduler runs physics and chemistry in a fixed order every tick (declared
//! by their layers: Physics before Chemistry). Within chemistry, cells are
//! processed in the sorted order of the active set. Nothing here reads a clock, a
//! wall time, or a RAM state; given the same world it does the same thing. The
//! reactor itself is bit-for-bit reproducible (proven in au-chem). So chemistry is
//! as deterministic as everything beneath it — a resumed world reacts identically
//! to one that never stopped.
//!
//! # Conservation
//!
//! Atoms are conserved by the reactor's structure (every event is atom-balanced,
//! populations never go negative). Energy is conserved because the heat the
//! reactor moves out of bonds is added, to the microjoule, into the shared thermal
//! field — and the ledger already accounts for every joule that crosses the
//! world's boundary via radiation and flux. Chemistry moves energy *within* the
//! world (bonds ⇄ thermal); it never crosses the boundary, so it touches no
//! ledger term. The invariant `E(now) − E(start) == in − out` is untouched by a
//! reaction, exactly as it should be: burning fuel does not change a sealed box's
//! total energy.

use au_core::event::kind;
use au_core::layer::Layer;
use au_core::schedule::{System, SystemDesc};
use au_data::chunk::ChunkCoord;
use au_physics::{
    derive, diffuse_explicit_perm, diffuse_implicit_perm, sweeps_for_deep_time, DiffuseScratch,
    Energy,
    GridSpec, MaterialRegistry, Phase,
};

use au_chem::{
    enumerate_reactions, react, CellChemistry, NetworkRules, ReactingCell, Reaction, SpeciesId,
};

use crate::chem_columns::{has_chemistry, species_column};
use crate::chemistry::Chemistry;
use crate::physics_columns::{has_physics, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};
use crate::world::World;

pub struct ChemistrySystem {
    chem: Chemistry,
    spec: GridSpec,
    n_species: usize,

    /// Open-ended mode (Phase 5a): reactions are generated, molecules are
    /// discovered. `None` runs the declared Phase 3 chemistry, byte-identically.
    open: Option<NetworkRules>,
    /// The reactions currently derivable from what the world contains. Rebuilt
    /// whenever the present-species set changes — a cache, not state: it is a
    /// pure function of the world, which is why resume needs no help from it.
    open_reactions: Vec<Reaction>,
    /// The species set the cache was built for.
    last_present: Vec<u32>,
    enumerated: bool,

    /// Derived catalysis (Phase 5c); `None` = off, and the world is 5b-identical.
    catalysis: Option<au_chem::CatalysisRules>,
    /// Ceiling on catalysed variants — a declared edge, see the config docs.
    max_catalysed: usize,

    /// Self-assembly; `None` = off, and the world is 5e-identical.
    membrane: Option<au_chem::MembraneRules>,
    /// Per-cell permeability, rebuilt each tick from the populations. Derived,
    /// never stored — there is no membrane column and nothing to serialize.
    perm: Vec<f64>,
    /// Amphiphilicity per species, extended as the world invents molecules. A
    /// molecule cannot change once interned, so a score is computed once and is
    /// right forever.
    species_amph: Vec<f64>,
    /// Scratch: the assembled (non-diffusing) part of one species, per cell.
    seq: Vec<i128>,

    /// Species transport (Phase 5b). D₀ of a 1-Da reference particle; 0 = off.
    diffusion_d0: f64,
    /// Deep-time switch for the diffusion solve — deliberately the SAME flag
    /// as conduction's (`physics.implicit_conduction`): one switch takes a
    /// whole world into deep time, heat and matter together.
    diff_implicit: bool,
    diff_iters: u32,
    diff_omega: f64,
    diff_scratch: DiffuseScratch,
    /// Per-cell mobility mask, rebuilt per chunk per tick: solutes move only
    /// where the host is fluid.
    mobile: Vec<bool>,
    /// Cached per-species rates D_s (m²/s), Graham-scaled from molecular mass.
    /// Grows with the registry; refreshed each tick before the chunk loop.
    species_d: Vec<f64>,

    // Working memory — all derived, none of it state, exactly as in physics.
    // Species populations, one scratch vec per species; plus the shared physics
    // fields we need to read (mass, material) and read-write (energy).
    pops: Vec<Vec<i128>>,
    mass: Vec<i128>,
    material: Vec<u16>,
    energy: Vec<i128>,

    /// Open boundaries (Phase 5h). Empty = a sealed world, byte-identical to 5g.
    ports: Vec<(usize, au_chem::Port)>,
    /// Atoms per element per species, rebuilt when the registry grows. Booking a
    /// port's delta means walking a molecular graph per species, and doing that
    /// every tick for every port would cost more than the chemistry.
    species_formula: Vec<Vec<i128>>,
    /// Scratch: one port's change, per species.
    port_delta: Vec<i128>,
    /// Scratch: this tick's boundary crossings, per species, summed over every
    /// port and chunk. Folded into the world's atom ledger once, after the chunk
    /// loop — the same discipline physics uses for energy, and for the same
    /// reason: a borrow of `world` inside the loop would forbid it.
    ///
    /// **Two vectors, not one signed vector.** A source and a sink acting on the
    /// same species cancel exactly at steady state, and a single accumulator
    /// would report a sealed world in the middle of a flow. See
    /// `AtomLedger::book_gross`.
    gained: Vec<i128>,
    lost: Vec<i128>,

    // Telemetry.
    pub last_events: u64,
    pub last_heat: Energy,
}

impl ChemistrySystem {
    pub fn from_config(world: &World) -> Option<(SystemDesc, ChemistrySystem)> {
        // The vocabulary is parsed once here; a malformed chemistry already
        // panicked in `World::new`, so by the time we are here it is either valid
        // or absent.
        let chem = match Chemistry::from_kv(world.config.as_map()) {
            Ok(Some(c)) => c,
            Ok(None) => return None,
            Err(e) => panic!("chemistry config invalid: {}", e),
        };
        let c = &world.config;
        let spec = GridSpec::new(
            c.u64_or("physics.grid.nx", 1) as usize,
            c.u64_or("physics.grid.ny", 1) as usize,
            c.u64_or("physics.grid.nz", 1) as usize,
            c.f64_or("physics.grid.cell_m", 10.0),
        );
        let n_species = chem.column_count();
        let open = chem.open;
        let chem_d0 = chem.diffusion_m2_s;
        let chem_membrane = chem.membrane;
        let chem_catalysis = chem.catalysis;
        let chem_max_catalysed = chem.max_catalysed;
        let chem_ports = chem.ports.clone();

        // Chemistry reads the shared thermal field that physics owns, so it must
        // run *after* physics each tick — declared by depending on the Physics
        // layer. The layer stack asserts Chemistry sits above Physics, so this is
        // legal (reading downward is allowed; only reading upward is refused).
        let desc = SystemDesc::new("chemistry.react", Layer::Chemistry)
            .reads(&[Layer::Universe, Layer::Physics, Layer::Chemistry])
            .every(1);

        Some((
            desc,
            ChemistrySystem {
                chem,
                spec,
                n_species,
                open,
                open_reactions: Vec::new(),
                last_present: Vec::new(),
                enumerated: false,
                catalysis: chem_catalysis,
                max_catalysed: chem_max_catalysed,
                membrane: chem_membrane,
                perm: Vec::new(),
                species_amph: Vec::new(),
                seq: Vec::new(),
                diffusion_d0: chem_d0,
                diff_implicit: world.config.bool_or("physics.implicit_conduction", false),
                diff_iters: world.config.u64_or("physics.conduction_iters", 24) as u32,
                diff_omega: world.config.f64_or("physics.conduction_omega", 1.0),
                diff_scratch: DiffuseScratch::new(),
                mobile: Vec::new(),
                species_d: Vec::new(),
                pops: vec![Vec::new(); n_species],
                mass: Vec::new(),
                material: Vec::new(),
                energy: Vec::new(),
                ports: chem_ports,
                species_formula: Vec::new(),
                port_delta: Vec::new(),
                gained: Vec::new(),
                lost: Vec::new(),
                last_events: 0,
                last_heat: Energy::ZERO,
            },
        ))
    }

    pub fn species_count(&self) -> usize {
        self.n_species
    }
    pub fn chemistry(&self) -> &Chemistry {
        &self.chem
    }
}

impl System<World> for ChemistrySystem {
    fn run(&mut self, world: &mut World) {
        if world.active.is_empty() || self.n_species == 0 {
            return;
        }
        let dt = world.clock.scale().as_secs_f64();
        if dt <= 0.0 {
            return;
        }

        let coords: Vec<ChunkCoord> = world.active.iter().copied().collect();
        for c in &coords {
            world.chunk(*c);
        }

        // ── Open chemistry: derive the reaction set from what exists. ───────
        //
        // "Present" is every species with a nonzero population anywhere in the
        // active set, plus the declared seeds (a seed at zero population adds
        // reactions of zero rate — harmless, and it keeps the set stable). The
        // cache is rebuilt only when that set changes; the enumeration itself
        // interns any product molecule the registry has never met, and each
        // interning is *history*, logged as a discovery event and thereby folded
        // into the world hash.
        if let Some(rules) = self.open {
            let mut present: std::collections::BTreeSet<u32> =
                (0..self.chem.species_count() as u32).collect();
            for c in &coords {
                let Some(chunk) = world.chunks.get(*c) else { continue };
                if !has_chemistry(chunk.columns(), self.n_species) {
                    continue;
                }
                for sp in 0..self.n_species {
                    let col = chunk.columns().get::<i128>(species_column(sp)).unwrap();
                    if col.as_slice().iter().any(|&v| v != 0) {
                        present.insert(sp as u32);
                    }
                }
            }
            let present: Vec<u32> = present.into_iter().collect();
            if !self.enumerated || present != self.last_present {
                let reg = &mut world.chem_registry;
                let before = reg.len();
                let ids: Vec<SpeciesId> = present.iter().map(|&i| SpeciesId(i)).collect();
                let rx = enumerate_reactions(&ids, reg, &self.chem.table, &self.chem.bond_model, &rules);
                let after = reg.len();
                // Shelf space: a reaction touching a species past the pool can
                // never be stocked, so it is dropped. The molecule stays in the
                // registry — the world *knows* it — it just cannot make it. A
                // declared edge, stated rather than hit.
                self.open_reactions = rx
                    .into_iter()
                    .filter(|r| {
                        r.products.iter().chain(r.reactants.iter())
                            .all(|t| (t.species.0 as usize) < self.n_species)
                    })
                    .collect();

                // ── Catalysis (Phase 5c). ───────────────────────────────────
                // The uncatalysed roads stay open — a catalyst makes a path
                // faster, it does not close the old one. What gets appended is
                // the extra paths, each with the catalyst on both sides, so its
                // rate is proportional to the catalyst's own concentration.
                // Where a catalyst happens to be the product of its own parent
                // reaction, the result reads `A + B + AB → 2 AB`. Nothing here
                // special-cases that; it is what the general rule says when the
                // structure lines up.
                if let Some(crules) = self.catalysis {
                    let base = self.open_reactions.clone();
                    let mut vars = au_chem::catalysed_variants(
                        &base,
                        &ids,
                        &world.chem_registry,
                        &self.chem.table,
                        &self.chem.bond_model,
                        &crules,
                        &rules,
                    );
                    if vars.len() > self.max_catalysed {
                        // Keep the strongest. Sorting by the achieved reduction
                        // (largest first) and breaking ties by ids keeps the
                        // truncation a pure function of the network.
                        vars.sort_by(|a, b| {
                            let ra = base[a.parent].activation_j - a.reaction.activation_j;
                            let rb = base[b.parent].activation_j - b.reaction.activation_j;
                            rb.partial_cmp(&ra)
                                .unwrap_or(std::cmp::Ordering::Equal)
                                .then(a.catalyst.0.cmp(&b.catalyst.0))
                                .then(a.parent.cmp(&b.parent))
                        });
                        vars.truncate(self.max_catalysed);
                        vars.sort_by(|a, b| {
                            a.parent.cmp(&b.parent).then(a.catalyst.0.cmp(&b.catalyst.0))
                        });
                    }
                    self.open_reactions.extend(vars.into_iter().map(|v| v.reaction));
                }
                for id in before..after {
                    let atoms =
                        world.chem_registry.get(SpeciesId(id as u32)).unwrap().atom_count();
                    world.emit(
                        Layer::Chemistry,
                        kind::CHEMISTRY_DISCOVERY,
                        id as u64,
                        atoms as u64,
                        after as u64,
                    );
                }
                self.last_present = present;
                self.enumerated = true;
            }
        }

        // ── Species mobilities, Graham-scaled from the discovered graphs. ───
        // D_s = D₀·√(1 Da / m_s): the reference scale is declared, the per-
        // species *ratios* are read off molecular structure — so a molecule the
        // world invented yesterday already knows how fast it moves today.
        if self.diffusion_d0 > 0.0 {
            self.species_d.resize(self.n_species, 0.0);
            for s_id in 0..self.n_species.min(world.chem_registry.len()) {
                let m = world.chem_registry.get(au_chem::SpeciesId(s_id as u32)).unwrap();
                let mass_mda = m.mass_mda(&self.chem.table).max(1) as f64;
                self.species_d[s_id] = self.diffusion_d0 * (1000.0 / mass_mda).sqrt();
            }
        }

        // ── Amphiphilicity, extended as the vocabulary grows. ───────────────
        // A molecule cannot change once interned, so a score computed today is
        // still right forever; only new ids need work.
        if self.membrane.is_some() {
            let known = world.chem_registry.len().min(self.n_species);
            while self.species_amph.len() < known {
                let id = au_chem::SpeciesId(self.species_amph.len() as u32);
                let m = world.chem_registry.get(id).unwrap();
                self.species_amph.push(au_chem::amphiphilicity(m, &self.chem.table));
            }
        }

        // ── Atoms per species, extended as the vocabulary grows. ────────────
        // Same argument as amphiphilicity: a molecule's formula is fixed the
        // moment it is interned, so only new ids need work. Needed only where
        // there are ports — a sealed world never books anything.
        if !self.ports.is_empty() {
            let known = world.chem_registry.len().min(self.n_species);
            while self.species_formula.len() < known {
                let id = au_chem::SpeciesId(self.species_formula.len() as u32);
                let m = world.chem_registry.get(id).unwrap();
                self.species_formula.push(m.formula(&self.chem.table));
            }
        }

        let mats: &MaterialRegistry = &world.materials;
        let chunks = &mut world.chunks;
        let reactions: &[Reaction] =
            if self.open.is_some() { &self.open_reactions } else { &self.chem.reactions };

        let mut total_events: u64 = 0;
        let mut net_heat: i128 = 0;
        self.gained.clear();
        self.gained.resize(self.n_species, 0);
        self.lost.clear();
        self.lost.resize(self.n_species, 0);

        for c in &coords {
            let Some(chunk) = chunks.get_mut(*c) else { continue };
            // Chemistry needs both its own species columns and the physics fields
            // it shares (mass, material, energy). A chunk lacking either is not a
            // chemistry chunk.
            if !has_physics(chunk.columns()) || !has_chemistry(chunk.columns(), self.n_species) {
                continue;
            }

            // ── Gather. ─────────────────────────────────────────────────────
            {
                let cols = chunk.columns();
                self.mass.clear();
                self.mass.extend_from_slice(cols.get::<i128>(PHYS_MASS).unwrap().as_slice());
                self.material.clear();
                self.material
                    .extend_from_slice(cols.get::<u16>(PHYS_MATERIAL).unwrap().as_slice());
                self.energy.clear();
                self.energy.extend_from_slice(cols.get::<i128>(PHYS_ENERGY).unwrap().as_slice());
                for s in 0..self.n_species {
                    self.pops[s].clear();
                    self.pops[s].extend_from_slice(
                        cols.get::<i128>(species_column(s)).unwrap().as_slice(),
                    );
                }
            }
            let cells = self.spec.cells();
            if self.mass.len() != cells {
                continue; // built for another grid; refuse to integrate garbage
            }

            // ── Compute, cell by cell. ──────────────────────────────────────
            // Each cell is an independent reactor: local populations, local
            // temperature (derived from the shared energy), local heat back into
            // the shared energy. No transport here — that is physics' job, and it
            // will move the heat and (eventually) the molecules on its own ticks.
            let volume = self.chem.cell_volume_m3;
            for cell in 0..cells {
                // Skip vacuum: no mass means no meaningful temperature and nothing
                // to react. (A cell can hold molecules only where there is matter.)
                let mass_kg = self.mass[cell] as f64 / au_physics::AG_PER_KG as f64;
                if mass_kg <= 0.0 {
                    continue;
                }
                // Any molecules here at all?
                let mut any = false;
                for s in 0..self.n_species {
                    if self.pops[s][cell] != 0 {
                        any = true;
                        break;
                    }
                }
                if !any {
                    continue;
                }

                // Temperature: derived from the shared energy, the SAME way physics
                // derives it, so the two never disagree. Material gives the heat
                // capacity; pressure is left at zero here (Phase 3b cells are the
                // gas-phase reactors of the demo; a full pressure coupling is a
                // Phase 4 refinement).
                let mat = mats
                    .get(au_physics::MaterialId(self.material[cell]))
                    .unwrap_or_else(|| mats.get(au_physics::MaterialId::VACUUM).unwrap());
                let energy_j = Energy(self.energy[cell]).as_joules();
                let state = derive(mass_kg, energy_j, mat, 0.0, volume);
                let temp = state.temperature;
                if temp <= 0.0 {
                    continue;
                }

                // Assemble a per-cell chemistry view over the scratch populations.
                // CellChemistry wants a contiguous population vector; we build one
                // for this cell, react it, and write the results back into the
                // per-species scratch.
                let mut cc = CellChemistry::new(self.n_species);
                for s in 0..self.n_species {
                    cc.set(au_chem::SpeciesId(s as u32), self.pops[s][cell]);
                }
                let mut thermal = Energy(self.energy[cell]);
                let report = {
                    let mut rc = ReactingCell {
                        chem: &mut cc,
                        thermal_energy: &mut thermal,
                        temperature_k: temp,
                        volume_m3: volume,
                    };
                    react(&mut rc, reactions, dt)
                };

                // Write back: populations into scratch, the (possibly changed)
                // energy into the shared thermal field. This is where the heat
                // enters the field physics conducts.
                for s in 0..self.n_species {
                    self.pops[s][cell] = cc.get(au_chem::SpeciesId(s as u32));
                }
                self.energy[cell] = thermal.0;

                total_events += report.events;
                net_heat += report.heat_released.0;
            }

            // ── Diffuse (Phase 5b): the discovered move. ────────────────────
            // Solutes flow down their own count gradients wherever the host is
            // fluid — the phase read from the same derive() the rates trust, so
            // "is this cell an ocean or an ice sheet" has exactly one answer.
            // A freezing cell exits the mask and traps its cargo where it sat.
            if self.diffusion_d0 > 0.0 && cells > 1 {
                let volume = self.chem.cell_volume_m3;
                self.mobile.clear();
                self.mobile.resize(cells, false);
                let mut ndim = 0;
                if self.spec.nx > 1 {
                    ndim += 1;
                }
                if self.spec.ny > 1 {
                    ndim += 1;
                }
                if self.spec.nz > 1 {
                    ndim += 1;
                }
                for cell in 0..cells {
                    let mass_kg = self.mass[cell] as f64 / au_physics::AG_PER_KG as f64;
                    if mass_kg <= 0.0 {
                        continue;
                    }
                    let mat = mats
                        .get(au_physics::MaterialId(self.material[cell]))
                        .unwrap_or_else(|| mats.get(au_physics::MaterialId::VACUUM).unwrap());
                    let st = derive(mass_kg, Energy(self.energy[cell]).as_joules(), mat, 0.0, volume);
                    self.mobile[cell] =
                        matches!(st.phase, Phase::Liquid | Phase::Gas);
                }

                // ── Membranes. ──────────────────────────────────────────────
                // Amphiphiles past their critical concentration are surface, not
                // solution, and surface is what the rest of the world has to get
                // through. A cell whose own chemistry builds surfactant thereby
                // builds its own barrier — nothing here decides that a cell
                // "should" become a compartment.
                let perm_opt: Option<&[f64]> = if let Some(mrules) = self.membrane {
                    self.perm.clear();
                    self.perm.resize(cells, 1.0);
                    for cell in 0..cells {
                        let mut assembled = 0i128;
                        for s in 0..self.n_species {
                            let score = self.species_amph.get(s).copied().unwrap_or(0.0);
                            assembled += au_chem::assembled_count(
                                self.pops[s][cell], score, volume, &mrules,
                            );
                        }
                        let cov = au_chem::coverage(assembled, volume, &mrules);
                        self.perm[cell] = au_chem::permeability(cov, &mrules);
                    }
                    Some(&self.perm[..])
                } else {
                    None
                };

                let dx = self.spec.cell_m;
                for s_id in 0..self.n_species {
                    let d = *self.species_d.get(s_id).unwrap_or(&0.0);
                    if d <= 0.0 || self.pops[s_id].iter().all(|&v| v == 0) {
                        continue;
                    }
                    let k = d * dt / (dx * dx);

                    // ── Assembly is a sink. ─────────────────────────────────
                    // A molecule in a membrane is part of a *structure*, not a
                    // solute: it has no independent gradient to run down. Only
                    // the free fraction diffuses, so it is subtracted here and
                    // added back afterwards — the same integers both ways, so
                    // conservation and non-negativity are untouched. Without
                    // this a compartment would dissolve itself by diffusion the
                    // moment it formed.
                    let sequestered = if let Some(mrules) = self.membrane {
                        let score = self.species_amph.get(s_id).copied().unwrap_or(0.0);
                        if score >= mrules.min_score {
                            self.seq.clear();
                            self.seq.resize(cells, 0);
                            for cell in 0..cells {
                                let a = au_chem::assembled_count(
                                    self.pops[s_id][cell], score, volume, &mrules,
                                );
                                self.seq[cell] = a;
                                self.pops[s_id][cell] -= a;
                            }
                            true
                        } else { false }
                    } else { false };
                    if self.diff_implicit {
                        // Sweeps scale with the grid span squared — the price
                        // of single-grid Gauss–Seidel at deep-time rates, paid
                        // knowingly (see au_physics::diffuse's module docs).
                        let sweeps = sweeps_for_deep_time(self.spec, self.diff_iters);
                        diffuse_implicit_perm(
                            self.spec,
                            &mut self.pops[s_id],
                            &self.mobile,
                            perm_opt,
                            k,
                            sweeps,
                            self.diff_omega,
                            &mut self.diff_scratch,
                        );
                    } else {
                        // Explicit path: substep to the stability ceiling, and
                        // refuse pathological cases loudly — a world that needs
                        // ten thousand diffusion substeps per tick has outgrown
                        // the explicit scheme and should say
                        // physics.implicit_conduction = true.
                        let k_max = 0.4 / (2.0 * ndim.max(1) as f64);
                        let n_sub = (k / k_max).ceil().max(1.0);
                        assert!(
                            n_sub <= 10_000.0,
                            "species diffusion needs {} explicit substeps per tick; \
                             enable physics.implicit_conduction for deep time",
                            n_sub
                        );
                        let h = k / n_sub;
                        for _ in 0..n_sub as u32 {
                            diffuse_explicit_perm(
                                self.spec,
                                &mut self.pops[s_id],
                                &self.mobile,
                                perm_opt,
                                h,
                                &mut self.diff_scratch,
                            );
                        }
                    }

                    // Put the membrane back where it was: the same integers that
                    // were removed, so the round trip is exact. (A previous
                    // attempt at this subtracted without restoring — it compiled
                    // perfectly and quietly destroyed matter on every tick.)
                    if sequestered {
                        for cell in 0..cells {
                            self.pops[s_id][cell] += self.seq[cell];
                        }
                    }
                }
            }

            // ── The boundary (Phase 5h). ────────────────────────────────────
            //
            // Last, and deliberately so. Reactions decide what the cell makes;
            // transport decides where it goes; the port then states what is true
            // at the boundary regardless of either. Running it first would let a
            // tick of chemistry consume the feed before anything else saw it, and
            // would leave a sink cell holding whatever diffused in after it was
            // emptied — an outflow that lags its own boundary condition by a tick.
            //
            // Each port is an exact integer clamp, so this is idempotent: the
            // second application in a row is a no-op. Every molecule it adds or
            // removes is booked, per element, into the world's atom ledger.
            if !self.ports.is_empty() {
                let cells = self.mass.len();
                for (cell, port) in &self.ports {
                    let cell = *cell;
                    if cell >= cells {
                        continue;
                    }
                    au_chem::port_delta(
                        port,
                        self.n_species,
                        |s| self.pops[s][cell],
                        &mut self.port_delta,
                    );
                    for (s, &d) in self.port_delta.iter().enumerate() {
                        if d > 0 {
                            self.pops[s][cell] += d;
                            self.gained[s] += d;
                        } else if d < 0 {
                            self.pops[s][cell] += d;
                            self.lost[s] += -d;
                        }
                    }
                }
            }

            // ── Scatter. `columns_mut` dirties the chunk — correct: a chunk
            // chemistry has touched is no longer a function of the seed and can
            // never be evicted again, exactly as with physics.
            {
                let cols = chunk.columns_mut();
                cols.get_mut::<i128>(PHYS_ENERGY)
                    .unwrap()
                    .as_mut_slice()
                    .copy_from_slice(&self.energy);
                for s in 0..self.n_species {
                    cols.get_mut::<i128>(species_column(s))
                        .unwrap()
                        .as_mut_slice()
                        .copy_from_slice(&self.pops[s]);
                }
            }
        }

        // Fold this tick's boundary crossings into the books. Once, after the
        // chunk loop, for the same reason physics does it that way: `world` is
        // borrowed inside the loop and cannot also be written there.
        if !self.ports.is_empty() {
            world.atoms.book_gross(&self.gained, &self.lost, &self.species_formula);
        }

        self.last_events = total_events;
        self.last_heat = Energy(net_heat);

        // Report into the causal history — but only if something happened, so a
        // dormant chemistry adds nothing to the hash. The event feeds the world
        // hash, which is exactly why a reaction must be deterministic: the same
        // universe must produce the same events.
        if total_events > 0 {
            world.emit(
                Layer::Chemistry,
                kind::CHEMISTRY_HEAT,
                net_heat.unsigned_abs() as u64,
                (net_heat > 0) as u64,
                total_events,
            );
        }
    }
}
