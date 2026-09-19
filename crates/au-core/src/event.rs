//! The causal record.
//!
//! VISION's Fundamental Rule is that everything exists because something caused
//! it. That claim is only worth anything if we can *check* it — if, when a
//! lineage goes extinct in year 12,000,000, we can walk backwards and find the
//! climate shift that did it.
//!
//! So significant state changes emit an event. Events are:
//!
//! * **Fixed-size and numeric.** 32 bytes. No strings, no allocation. At the
//!   scales this project is aiming for, a `String` per event would eat the heap
//!   long before anything interesting evolved.
//! * **Append-only.** History is not editable.
//! * **Sparse.** ARCHITECTURE §19: store important events, regenerate the rest.
//!   A tick where nothing notable happens writes nothing at all.
//!
//! `kind` is a per-layer opcode. The interpretation of `a`/`b`/`c` is defined
//! by the layer that emits it — the core neither knows nor cares.

use crate::hash::{Hasher, WorldHash};
use crate::layer::Layer;
use crate::time::Tick;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Event {
    pub tick: Tick,
    pub layer: Layer,
    pub kind: u16,
    pub a: u64,
    pub b: u64,
    pub c: u64,
}

/// Event kinds owned by the Universe layer. Each layer gets its own module of
/// these as it is built.
pub mod kind {
    /// The simulation clock's seconds-per-tick changed. a = new secs, b = new frac.
    pub const TIME_SCALE_CHANGED: u16 = 0x0001;

    /// The physics could not resolve the timestep the clock asked for, and the
    /// clock has been slowed to what diffusion actually permits.
    /// a = capped secs/tick, b = capped frac, c = the substep budget it burned.
    ///
    /// **This is the engine admitting a limit, not reporting an error.** You
    /// cannot fast-forward a diffusion equation; deep time needs coarse-grained
    /// physics, not the same physics run recklessly.
    pub const PHYSICS_TIME_CAPPED: u16 = 0x0100;

    /// A chemistry step released (or absorbed) a net quantity of heat into the
    /// shared thermal field of the active set.
    /// a = |net heat| in microjoules, b = 1 if the cell(s) were heated (net
    /// exothermic) else 0, c = total reaction events fired this step.
    ///
    /// This is chemistry reporting into the same causal history physics does —
    /// the coupling made visible. It feeds the world hash, so a reaction that
    /// fires in one run must fire identically in every run.
    pub const CHEMISTRY_HEAT: u16 = 0x0200;

    /// An open chemistry interned a molecule the world had never seen.
    /// a = the new species id, b = its atom count, c = total species known now.
    ///
    /// Discovery is *history* — the moment a substance began to exist for this
    /// world — so it goes in the log like any other cause, and it feeds the
    /// world hash: two worlds that found different molecules are different
    /// worlds even before any population differs.
    pub const CHEMISTRY_DISCOVERY: u16 = 0x0201;

    // NOTE, deliberately: there is no CHUNK_EVICTED or CHUNK_GENERATED here.
    //
    // Memory management is not causality. The event log feeds the world hash,
    // so an eviction event would make *how much RAM the host had* part of the
    // identity of the universe — two machines would compute different
    // histories. Chunk residency is telemetry (see `ChunkStore::stat_*`) and
    // must remain invisible to the simulation.
}

impl WorldHash for Event {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u64(self.tick.0);
        h.write_u8(self.layer as u8);
        h.write_u16(self.kind);
        h.write_u64(self.a);
        h.write_u64(self.b);
        h.write_u64(self.c);
    }
}

/// Append-only history.
///
/// `cap` bounds live memory: once full, the oldest events are dropped from RAM.
/// This is not data loss — the log is meant to be streamed to disk behind the
/// simulation. The cap exists so that a world left running for a week does not
/// die of its own diary.
#[derive(Clone, Debug)]
pub struct EventLog {
    events: Vec<Event>,
    cap: usize,
    /// Total ever emitted, including those aged out. Part of world identity, so
    /// that a truncated log is still a *detectably* consistent one.
    total: u64,
}

impl EventLog {
    pub fn new(cap: usize) -> Self {
        EventLog { events: Vec::new(), cap, total: 0 }
    }

    pub fn emit(&mut self, e: Event) {
        if self.cap > 0 && self.events.len() == self.cap {
            self.events.remove(0);
        }
        self.events.push(e);
        self.total += 1;
    }

    pub fn events(&self) -> &[Event] {
        &self.events
    }

    pub fn total(&self) -> u64 {
        self.total
    }

    /// The in-memory bound. Must be persisted: a log restored with the wrong
    /// capacity silently starts discarding history.
    pub fn cap(&self) -> usize {
        self.cap
    }

    pub fn len(&self) -> usize {
        self.events.len()
    }

    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }

    /// Every event a given layer emitted. The first tool of the debugger.
    pub fn by_layer(&self, layer: Layer) -> impl Iterator<Item = &Event> {
        self.events.iter().filter(move |e| e.layer == layer)
    }

    pub fn restore(events: Vec<Event>, cap: usize, total: u64) -> Self {
        EventLog { events, cap, total }
    }
}

impl WorldHash for EventLog {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u64(self.total);
        h.write_u64(self.events.len() as u64);
        for e in &self.events {
            e.hash_into(h);
        }
    }
}
