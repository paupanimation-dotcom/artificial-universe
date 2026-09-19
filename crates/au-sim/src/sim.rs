//! The tick loop.
//!
//! Small on purpose. This file should still be about this size when there are
//! thirteen layers in it, because layers are added by *registering systems*,
//! not by editing the loop. If `tick()` ever grows a `match` on layer, or a
//! special case for organisms, the architecture has failed and this is where
//! it will show first.

use au_core::hash::WorldHash;
use au_core::layer::Layer;
use au_core::schedule::Scheduler;
use au_core::time::Tick;
use au_data::snapshot;

use crate::config::Config;
use au_core::ids::EntityAllocator;

use crate::systems::chemistry::ChemistrySystem;
use crate::systems::protocell::ProtocellSystem;
use crate::systems::physics::PhysicsSystem;
use crate::systems::time_scale::TimeScalePolicy;
use crate::world::World;

pub struct Simulation {
    pub world: World,
    scheduler: Scheduler<World>,
    evict_every: u64,
}


/// If the config declares a planet, install its generator on the world — so a
/// planet world's untouched chunks are the planet, automatically, in both a fresh
/// boot and a resume. (A test or demo can still override with its own generator
/// afterwards; last set wins, as before.)
fn install_planet(world: &mut World) {
    let n_species = match crate::chemistry::Chemistry::from_kv(world.config.as_map()) {
        Ok(Some(c)) => c.species_count(),
        _ => 0,
    };
    if let Some(gen) = crate::planet_gen::PlanetGenerator::from_config(
        world.seed,
        &world.config,
        &world.materials,
        n_species,
    ) {
        world.set_generator(Box::new(gen));
    }
}

/// Marks a save frame that carries an atom ledger. Read from the last eight
/// bytes; absent in every world where nothing crossed a boundary.
const ATOM_TRAILER: &[u8; 8] = b"AUATOMS1";

/// Marks a save frame carrying protocells and the entity allocator.
const LIFE_TRAILER: &[u8; 8] = b"AULIFE01";

impl Simulation {
    /// Build a universe from configuration.
    ///
    /// Registration order does not matter: the scheduler sorts by layer, so the
    /// dependency chain in ARCHITECTURE §1 is enforced regardless of how
    /// carelessly systems are added here.
    pub fn new(config: Config) -> Simulation {
        let evict_every = config.u64_or("chunks.evict_every_ticks", 0);
        let mut world = World::new(config);
        let mut scheduler = Scheduler::new();

        // --- Universe layer -------------------------------------------------
        if let Some((desc, sys)) = TimeScalePolicy::from_config(&world) {
            scheduler.must_register(desc, sys);
        }

        // --- Physics layer --------------------------------------------------
        // Thermodynamics. Also, quietly, the authority on how fast the clock is
        // allowed to run — because diffusion has a speed limit and the Universe
        // layer above does not know that.
        if let Some((desc, sys)) = PhysicsSystem::from_config(&world) {
            scheduler.must_register(desc, sys);
        }

        // --- Chemistry layer ------------------------------------------------
        // Reactions run here, reading the shared thermal field physics owns —
        // which is why this must sit above Physics in the layer stack. An
        // exothermic reaction warms the field physics then conducts and radiates.
        // Present only if the world declares a chemistry in config.
        if let Some((desc, sys)) = ChemistrySystem::from_config(&world) {
            scheduler.must_register(desc, sys);
        }

        // --- Life chemistry layer -------------------------------------------
        // Protocells: bags that nucleate out of assembled surfactant, feed
        // through their own membranes, run their own chemistry at their own
        // concentration, and divide when geometry allows. Present only if the
        // world asks for them.
        if let Some((desc, sys)) = ProtocellSystem::from_config(&world) {
            scheduler.must_register(desc, sys);
        }

        // --- Matter, Planet, ... --------------------------------------------
        // Later phases. Empty is the correct state; a stub that returns plausible
        // numbers is worse than nothing, because everything above it would be
        // built against a fiction and nobody would find out for a year.

        install_planet(&mut world);

        Simulation { world, scheduler, evict_every }
    }

    /// Advance the universe by exactly one tick.
    pub fn tick(&mut self) {
        // Systems run for the tick they can see. The clock advances afterwards,
        // so within a tick every system agrees on what time it is — otherwise
        // a system's behaviour would depend on where it sat in the schedule,
        // which is exactly the kind of hidden ordering dependency that breaks
        // determinism the day we go parallel.
        let t = self.world.clock.tick();
        self.scheduler.run_tick(t, &mut self.world);

        // Housekeeping, NOT a simulation event.
        //
        // I wrote this as `world.emit(CHUNK_EVICTED, ...)` first, and it was a
        // real bug: the event log feeds the world hash, so emitting here would
        // have made *memory pressure* part of the identity of the universe. Two
        // machines with different RAM would have computed different histories,
        // and the failure would have surfaced as "evolution is subtly different
        // on the server", which is a nightmare to trace.
        //
        // The invariant, stated plainly: **nothing the engine does to manage
        // memory may be observable to the simulation.** Eviction is counted in
        // telemetry and is invisible to the world hash. There is a test.
        if self.evict_every > 0 && t.0 % self.evict_every == 0 && t.0 > 0 {
            self.world.evict_clean();
        }

        self.world.clock.advance();
    }

    pub fn run(&mut self, ticks: u64) {
        for _ in 0..ticks {
            self.tick();
        }
    }

    /// Run until a given tick, doing nothing if already past it.
    pub fn run_until(&mut self, target: Tick) {
        while self.world.clock.tick() < target {
            self.tick();
        }
    }

    /// The identity of this universe right now. Two runs of the same seed must
    /// agree on this at every tick, forever.
    pub fn hash(&self) -> u64 {
        self.world.world_hash()
    }

    pub fn save(&self) -> Vec<u8> {
        let active: Vec<_> = self.world.active.iter().copied().collect();
        let l = self.world.ledger;
        let snap = snapshot::encode(
            self.world.seed,
            &self.world.clock,
            &self.world.chunks,
            &self.world.events,
            &active,
            [
                l.energy_in.0,
                l.energy_out.0,
                l.mass_in.0,
                l.mass_out.0,
                l.impulse[0].0,
                l.impulse[1].0,
                l.impulse[2].0,
                0,
            ],
        );
        // The save frame: [u64 snapshot length][snapshot][species registry].
        //
        // Under a declared chemistry the registry is a pure function of config
        // and never needed saving. Discovery changed its nature: which molecules
        // exist is now something this world *did*, and a resumed world must be
        // able to name every species its cells hold populations of. Framing the
        // registry here keeps au-data's codec untouched — the snapshot stays
        // exactly the format it was, wrapped rather than revised.
        let reg = self.world.chem_registry.to_bytes();
        let mut out = Vec::with_capacity(8 + snap.len() + reg.len());
        out.extend_from_slice(&(snap.len() as u64).to_le_bytes());
        out.extend_from_slice(&snap);
        out.extend_from_slice(&reg);

        // ── The atom trailer (Phase 5h). ────────────────────────────────────
        //
        // Once matter can cross the boundary, "how many atoms are there" stops
        // being a constant and becomes `start + in - out`. The in and out are
        // history, so they must survive a save — and per element, because a bug
        // that turned carbon into oxygen would balance in bulk mass and still be
        // a catastrophe.
        //
        // Appended as a *trailer*, read backwards from the end, rather than
        // inserted into the frame. The registry section runs to the end of the
        // frame and carries no length, so anything placed after it would have to
        // be told apart from registry bytes; anything placed before it would move
        // the registry and break every save in existence. A trailer identified by
        // a magic at the very end costs one comparison, leaves the frame layout
        // untouched, and lets a world written before this phase load with an
        // empty ledger — which is the correct reading of a save from a world
        // where nothing could cross.
        // Order matters and is the whole reason the reader works: the atom
        // trailer is written first and the life trailer last, so `load` peels
        // life, then atoms, then finds the registry. Written the other way round
        // (as the first draft was) the reader looks for a life magic that is not
        // at the end and quietly fails to parse a living world.
        if self.world.atoms.any_flow() {
            // Length last, immediately before the magic: the reader works
            // backwards from the end of the frame and cannot find the start of
            // the trailer until it knows how long the trailer is.
            for v in self.world.atoms.inflow.iter().chain(self.world.atoms.outflow.iter()) {
                out.extend_from_slice(&v.to_le_bytes());
            }
            out.extend_from_slice(&(self.world.atoms.len() as u64).to_le_bytes());
            out.extend_from_slice(ATOM_TRAILER);
        }
        // ── The life trailer (Phase 5i step 2). ─────────────────────────────
        //
        // Protocells are individuals, so their contents are state and their ids
        // are identity. The entity allocator goes with them: a resumed world with
        // a reset allocator would reissue ids that living bags already hold.
        //
        // Written before the atom trailer and read after it, so the frame peels
        // from the end one section at a time. Absent when there are no
        // protocells, which keeps every pre-5i save loadable.
        if !self.world.protocells.is_empty() || self.world.entities.live_count() > 0 {
            let body = self.world.protocells.to_bytes();
            let (gens, free) = self.world.entities.raw();
            out.extend_from_slice(&body);
            for g in gens.iter() {
                out.extend_from_slice(&g.to_le_bytes());
            }
            for f in free.iter() {
                out.extend_from_slice(&f.to_le_bytes());
            }
            // Every length last, so the whole section can be located walking
            // backwards from the magic. The first draft put the allocator's
            // lengths at the *front* of the section and then indexed from the
            // front of the remaining bytes — which is the registry, not the
            // store. It parsed garbage and refused to load a living world.
            out.extend_from_slice(&(body.len() as u64).to_le_bytes());
            out.extend_from_slice(&(gens.len() as u64).to_le_bytes());
            out.extend_from_slice(&(free.len() as u64).to_le_bytes());
            out.extend_from_slice(LIFE_TRAILER);
        }

        out
    }

    /// Resurrect a universe. Must produce a simulation bit-identical to one
    /// that never stopped — there is a test that proves it.
    pub fn load(bytes: &[u8], config: Config) -> Result<Simulation, au_data::CodecError> {
        let evict_every = config.u64_or("chunks.evict_every_ticks", 0);
        // Registry must be rebuilt before decode: a save names its columns by
        // number, and an unregistered number is an error, not a shrug.
        let mut probe = World::new(config.clone());
        register_columns(&mut probe);
        // Split the frame `save` wrote: [len][snapshot][registry].
        if bytes.len() < 8 {
            return Err(au_data::CodecError::UnexpectedEof { needed: 8, had: bytes.len() });
        }
        let snap_len = u64::from_le_bytes(bytes[..8].try_into().unwrap()) as usize;
        if bytes.len() < 8 + snap_len {
            return Err(au_data::CodecError::UnexpectedEof {
                needed: 8 + snap_len,
                had: bytes.len(),
            });
        }
        let mut reg_bytes = &bytes[8 + snap_len..];

        // Peel the life trailer first — it was written last. See `save`.
        let mut protocells = crate::protocell::ProtocellStore::new();
        let mut entities = EntityAllocator::new();
        if reg_bytes.len() >= LIFE_TRAILER.len()
            && &reg_bytes[reg_bytes.len() - LIFE_TRAILER.len()..] == LIFE_TRAILER
        {
            let body = &reg_bytes[..reg_bytes.len() - LIFE_TRAILER.len()];
            if body.len() < 24 {
                return Err(au_data::CodecError::UnexpectedEof { needed: 24, had: body.len() });
            }
            let n = body.len();
            let store_len = u64::from_le_bytes(body[n - 24..n - 16].try_into().unwrap()) as usize;
            let ng = u64::from_le_bytes(body[n - 16..n - 8].try_into().unwrap()) as usize;
            let nf = u64::from_le_bytes(body[n - 8..].try_into().unwrap()) as usize;
            let section = store_len + ng * 4 + nf * 4 + 24;
            if n < section {
                return Err(au_data::CodecError::Invalid("life trailer"));
            }
            let base = n - section;
            let gens: Vec<u32> = (0..ng)
                .map(|i| {
                    let o = base + store_len + i * 4;
                    u32::from_le_bytes(body[o..o + 4].try_into().unwrap())
                })
                .collect();
            let free: Vec<u32> = (0..nf)
                .map(|i| {
                    let o = base + store_len + ng * 4 + i * 4;
                    u32::from_le_bytes(body[o..o + 4].try_into().unwrap())
                })
                .collect();
            protocells =
                crate::protocell::ProtocellStore::from_bytes(&body[base..base + store_len])
                    .ok_or(au_data::CodecError::Invalid("protocell store"))?;
            entities = EntityAllocator::restore(gens, free);
            reg_bytes = &reg_bytes[..base];
        }

        // Split off the atom trailer, if this world had ports. See `save`.
        let mut atoms = au_chem::AtomLedger::default();
        if reg_bytes.len() >= ATOM_TRAILER.len()
            && &reg_bytes[reg_bytes.len() - ATOM_TRAILER.len()..] == ATOM_TRAILER
        {
            let body = &reg_bytes[..reg_bytes.len() - ATOM_TRAILER.len()];
            if body.len() < 8 {
                return Err(au_data::CodecError::UnexpectedEof { needed: 8, had: body.len() });
            }
            let n = u64::from_le_bytes(body[body.len() - 8..].try_into().unwrap()) as usize;
            let need = 8 + 2 * n * 16;
            if body.len() < need {
                return Err(au_data::CodecError::UnexpectedEof { needed: need, had: body.len() });
            }
            let start = body.len() - need;
            let mut off = start;
            let read = |off: &mut usize| {
                let v = i128::from_le_bytes(body[*off..*off + 16].try_into().unwrap());
                *off += 16;
                v
            };
            let inflow: Vec<i128> = (0..n).map(|_| read(&mut off)).collect();
            let outflow: Vec<i128> = (0..n).map(|_| read(&mut off)).collect();
            atoms = au_chem::AtomLedger { inflow, outflow };
            reg_bytes = &reg_bytes[..start];
        }
        let snap = snapshot::decode(&bytes[8..8 + snap_len], &probe.registry)?;

        let active: std::collections::BTreeSet<_> = snap.active.iter().copied().collect();
        let c = snap.conserved;
        let ledger = au_physics::Ledger {
            energy_in: au_physics::Energy(c[0]),
            energy_out: au_physics::Energy(c[1]),
            mass_in: au_physics::Mass(c[2]),
            mass_out: au_physics::Mass(c[3]),
            impulse: [
                au_physics::Momentum(c[4]),
                au_physics::Momentum(c[5]),
                au_physics::Momentum(c[6]),
            ],
        };
        let (seed, clock, chunks, events) = snap.into_parts();
        let registry = probe.registry.clone();
        let world =
            World::restore(
                config, seed, clock, chunks, events, registry, active, ledger, atoms,
                protocells, entities,
            );

        let mut sim = Simulation {
            world,
            scheduler: Scheduler::new(),
            evict_every,
        };
        // Rebuild the schedule exactly as `new()` would. Systems hold no history
        // — only scratch buffers — so a rebuilt schedule is an identical one.
        // There is a test that proves it.
        if let Some((desc, sys)) = TimeScalePolicy::from_config(&sim.world) {
            sim.scheduler.must_register(desc, sys);
        }
        if let Some((desc, sys)) = PhysicsSystem::from_config(&sim.world) {
            sim.scheduler.must_register(desc, sys);
        }
        if let Some((desc, sys)) = ProtocellSystem::from_config(&sim.world) {
            sim.scheduler.must_register(desc, sys);
        }
        if let Some((desc, sys)) = ChemistrySystem::from_config(&sim.world) {
            sim.scheduler.must_register(desc, sys);
        }
        install_planet(&mut sim.world);
        // The discovered species come back by name, not by re-derivation —
        // discovery is history, and history is loaded, never replayed.
        sim.world.chem_registry = au_chem::SpeciesRegistry::from_bytes(reg_bytes)
            .map_err(|_| au_data::CodecError::Invalid("species registry"))?;
        Ok(sim)
    }

    pub fn schedule_report(&self) -> Vec<(&'static str, Layer, u64, u64)> {
        self.scheduler.report()
    }
}

/// Every column every layer owns.
///
/// Registration now happens in `World::new`, so this is a no-op kept for the
/// load path's readability. Physics owns `0x1xxx`; chemistry will own `0x2xxx`.
fn register_columns(_world: &mut World) {}
