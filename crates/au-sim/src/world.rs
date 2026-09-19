//! The World — the single mutable thing every system is handed.
//!
//! Note what is *not* here: planets, matter, organisms. Not because they are
//! coming later as fields on this struct, but because they will arrive as
//! **columns in chunks** and **systems in the scheduler**. `World` is the
//! substrate; it should stay roughly this size forever. If it starts growing a
//! field per layer, the architecture has quietly become a monolith and we
//! should stop and fix it.

use au_core::event::{Event, EventLog};
use au_core::hash::{Hasher, WorldHash};
use au_core::ids::EntityAllocator;
use au_core::layer::Layer;
use au_core::rng::{Domain, Rng};
use au_core::time::{SimClock, SimDuration, Tick};
use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, ChunkStore, EmptyGenerator};
use au_data::column::ColumnRegistry;
use au_physics::{Ledger, MaterialRegistry};
use std::collections::BTreeSet;

use crate::config::Config;

pub struct World {
    /// The number the whole universe unfolds from.
    pub seed: u64,
    pub clock: SimClock,
    pub chunks: ChunkStore,
    pub entities: EntityAllocator,
    pub events: EventLog,
    pub registry: ColumnRegistry,
    pub config: Config,

    /// The molecules this world knows. Seeded from config at boot; in an
    /// open chemistry it **grows** as reactions discover new graphs, which makes
    /// it history — serialized in every save, folded into the world hash.
    pub chem_registry: au_chem::SpeciesRegistry,

    /// The substances that exist. Empty on ship — Phase 3 will *derive* these
    /// from elements rather than list them. Config-driven so that a reloaded
    /// world runs on exactly the physics it was computed with.
    pub materials: MaterialRegistry,

    /// **The chunks under simulation.**
    ///
    /// Not "the chunks in memory". The distinction is the whole ballgame: if
    /// physics ran on whatever happened to be resident, then the set of things
    /// being simulated would be the set of things somebody *looked at*, and the
    /// universe would evolve differently depending on where the camera was
    /// pointed. Observation would become interaction.
    ///
    /// So this is explicit, snapshotted and hashed. [`World::activate`] is the
    /// only door in, and looking at a chunk does not knock on it.
    pub active: BTreeSet<ChunkCoord>,

    /// Every joule and gram that has crossed the boundary of the world.
    ///
    /// A planet is an open system — that is the only reason life is possible —
    /// so we cannot assert that the total never changes. We assert the stronger
    /// thing: `E(now) − E(start) == in − out`, exactly, as integers. See
    /// `au_physics::quantity`.
    pub ledger: Ledger,

    /// Every atom that has crossed the boundary of the world, per element.
    ///
    /// The twin of `ledger`, and separate from it for a reason: mass is a bulk
    /// quantity and would balance perfectly through a bug that turned carbon
    /// into oxygen. Chemistry needs the identity to hold *per element*, so this
    /// is counted per element. Empty until a world has ports (Phase 5h); a
    /// closed world carries a zero ledger forever and the old invariant —
    /// "totals never change" — is exactly the case `in - out == 0`.
    pub atoms: au_chem::AtomLedger,

    /// Every protocell in the world (Phase 5i step 2).
    ///
    /// **Species populations now live in two places.** The chunk columns hold
    /// the bulk medium; this holds what individuals have taken up. Any
    /// accounting of what the world contains must sum both — the first time in
    /// this project that "how much is there" has needed more than one place to
    /// look, and worth flagging because a check that forgets protocells will see
    /// matter vanish every time a vesicle buds.
    pub protocells: crate::protocell::ProtocellStore,

    generator: Box<dyn ChunkGenerator>,
}

impl World {
    pub fn new(config: Config) -> Self {
        let seed = config.u64_or("world.seed", 0);
        let spt = config.f64_or("clock.seconds_per_tick", 1.0);
        let cap = config.u64_or("events.capacity", 100_000) as usize;

        let materials = MaterialRegistry::from_kv(config.as_map())
            .unwrap_or_else(|e| panic!("bad material table: {}", e));

        let mut registry = ColumnRegistry::new();
        crate::physics_columns::register(&mut registry);
        // Chemistry columns, if the world declares a chemistry. The species count
        // is fixed at boot (see chem_columns), so the columns can be registered
        // now and round-trip through the snapshot codec like any other. A
        // malformed chemistry is a hard error here, not a silent skip — every
        // cell's populations would otherwise be meaningless.
        let (n_chem_cols, chem_registry, n_elements) =
            match crate::chemistry::Chemistry::from_kv(config.as_map()) {
                Ok(Some(c)) => (c.column_count(), c.species, c.table.len()),
                Ok(None) => (0, au_chem::SpeciesRegistry::new(), 0),
                Err(e) => panic!("bad chemistry table: {}", e),
            };
        crate::chem_columns::register(&mut registry, n_chem_cols);

        World {
            seed,
            clock: SimClock::new(SimDuration::from_secs_f64(spt)),
            chunks: ChunkStore::new(seed),
            entities: EntityAllocator::new(),
            events: EventLog::new(cap),
            registry,
            materials,
            chem_registry,
            active: BTreeSet::new(),
            ledger: Ledger::default(),
            atoms: au_chem::AtomLedger::new(n_elements),
            protocells: crate::protocell::ProtocellStore::new(),
            config,
            // Phase 1: the world is genuinely empty. Matter does not exist yet
            // (Phase 3). We are not going to sprinkle in placeholder terrain so
            // that the demo looks busier — that is exactly the "content" the
            // Golden Rule forbids, and it would be load-bearing before anyone
            // noticed.
            generator: Box::new(EmptyGenerator),
        }
    }

    /// Swap the world generator. Phase 4 (planets) installs the real one.
    ///
    /// Generators must be **pure** — `f(seed, coord)` and nothing else — or
    /// chunk eviction stops being lossless and the memory strategy dies.
    pub fn set_generator(&mut self, gen: Box<dyn ChunkGenerator>) {
        self.generator = gen;
    }

    /// Derive a private random stream. There is no global RNG; see `au_core::rng`.
    #[inline]
    pub fn rng(&self, domain: Domain, key: u64) -> Rng {
        // The tick is folded in so that the same system, at the same place, on
        // a different tick, gets a different stream — while still being a pure
        // function of (seed, where, when). Reproducible, and safe to run in
        // parallel because nothing is shared.
        Rng::derive(self.seed, domain, key, self.clock.tick().0)
    }

    #[inline]
    pub fn tick(&self) -> Tick {
        self.clock.tick()
    }

    pub fn emit(&mut self, layer: Layer, kind: u16, a: u64, b: u64, c: u64) {
        let tick = self.clock.tick();
        self.events.emit(Event { tick, layer, kind, a, b, c });
    }

    /// Touch a chunk, generating it if it has never been seen.
    ///
    /// **This does not enrol the chunk in the simulation.** Looking is not
    /// touching. See [`World::activate`].
    pub fn chunk(&mut self, coord: ChunkCoord) -> &mut Chunk {
        self.chunks.get_or_generate(coord, self.generator.as_ref())
    }

    /// Put a chunk under simulation.
    ///
    /// An explicit act, and a consequential one: from here on, physics runs in
    /// this chunk, which means it will be modified, which means it can never be
    /// evicted again. An active chunk costs memory forever.
    ///
    /// That is not a defect. It is the honest price of an evolving world, and it
    /// is precisely why level of detail is not optional: most of a planet must
    /// *not* be finely evolving, or the planet will not fit in the machine.
    ///
    /// In Phase 4 the planet decides what is active. Nothing else may.
    pub fn activate(&mut self, coord: ChunkCoord) {
        self.chunk(coord);
        self.active.insert(coord);
    }

    /// Drop everything that can be rebuilt from the seed.
    pub fn evict_clean(&mut self) -> usize {
        self.chunks.evict_clean()
    }

    #[allow(clippy::too_many_arguments)]
    pub fn restore(
        config: Config,
        seed: u64,
        clock: SimClock,
        chunks: ChunkStore,
        events: EventLog,
        registry: ColumnRegistry,
        active: BTreeSet<ChunkCoord>,
        ledger: Ledger,
        atoms: au_chem::AtomLedger,
        protocells: crate::protocell::ProtocellStore,
        entities: EntityAllocator,
    ) -> Self {
        let cfg_map = config.as_map();
        let materials = MaterialRegistry::from_kv(config.as_map())
            .unwrap_or_else(|e| panic!("bad material table: {}", e));
        // A save from a world where nothing crossed carries no atom trailer, so
        // `atoms` arrives empty. An empty ledger and a zeroed one are the same
        // statement about history and *must* be the same world: size it here, the
        // way `new` does, or a resumed sealed world hashes differently from the
        // one that never stopped purely because of how many elements it forgot to
        // write down. (Found by `catalysed_worlds_are_reproducible_and_resumable`,
        // which is exactly the test that should have found it.)
        let atoms = if atoms.is_empty() {
            let n = match crate::chemistry::Chemistry::from_kv(cfg_map) {
                Ok(Some(c)) => c.table.len(),
                _ => 0,
            };
            au_chem::AtomLedger::new(n)
        } else {
            atoms
        };
        World {
            seed,
            clock,
            chunks,
            // Restored, not fresh: protocells allocate ids, and a resumed world
            // that reset its allocator would hand out ids that already belong to
            // living bags — the exact stranger-in-the-family-tree failure the
            // generation counter was added in Phase 1 to prevent.
            entities,
            events,
            registry,
            materials,
            chem_registry: match crate::chemistry::Chemistry::from_kv(cfg_map) {
                Ok(Some(c)) => c.species,
                _ => au_chem::SpeciesRegistry::new(),
            },
            active,
            ledger,
            atoms,
            protocells,
            config,
            generator: Box::new(EmptyGenerator),
        }
    }
}

impl WorldHash for World {
    /// The identity of the universe, in 64 bits.
    ///
    /// Everything that can differ between two runs must be folded in here, and
    /// nothing that *cannot* differ should be — a hash that includes, say, the
    /// resident chunk count would report false divergence every time memory
    /// pressure differed. Note `ChunkStore` deliberately hashes only *dirty*
    /// chunks, for exactly this reason.
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u64(self.seed);
        h.write_u64(self.clock.tick().0);
        h.write_u64(self.clock.now().secs);
        h.write_u64(self.clock.now().frac);
        h.write_u64(self.clock.scale().secs);
        h.write_u64(self.clock.scale().frac);
        self.entities.hash_into(h);
        self.chunks.hash_into(h);
        self.events.hash_into(h);

        // The active set is world state, not a memory detail: it decides what is
        // simulated, and two worlds simulating different regions are different
        // worlds.
        h.write_u64(self.active.len() as u64);
        for c in &self.active {
            h.write_i32(c.x);
            h.write_i32(c.y);
            h.write_i32(c.z);
        }
        self.ledger.hash_into(h);
        self.atoms.hash_into(h);
        self.protocells.hash_into(h);

        // Discovered chemistry is history. Two worlds that have found different
        // molecules are different worlds even before a single population
        // differs — the registry's canonical forms, in id order, say so.
        h.write_u64(self.chem_registry.len() as u64);
        for (id, _) in self.chem_registry.iter() {
            h.write_bytes(&self.chem_registry.canonical(id).unwrap().0);
        }
    }
}
