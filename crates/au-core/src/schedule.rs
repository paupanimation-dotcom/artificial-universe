//! The multi-rate scheduler.
//!
//! ARCHITECTURE §3 asks for systems that update at wildly different rates —
//! physics constantly, chemistry frequently, evolution over generations,
//! civilizations over centuries. A single `update()` loop cannot express that,
//! and a loop that tries will spend 99.99% of its budget re-running ecology
//! for a world where nothing has changed.
//!
//! So systems don't "run each frame". They **register a period**, in ticks, and
//! the scheduler dispatches only what is due.
//!
//! This one mechanism buys three things at once:
//!
//! 1. **Multi-scale time.** Chemistry at period 1, evolution at period 65536.
//! 2. **Temporal LOD** (ARCHITECTURE §20). A distant region registers the same
//!    logical layer with a longer period and a coarser implementation. Nothing
//!    else in the codebase needs to know that "far away" exists.
//! 3. **Deterministic order.** Systems are sorted by (layer, registration
//!    index) — never by hash-map order, never by thread completion.
//!
//! `phase` staggers systems that share a period so they don't all land on the
//! same tick and produce a sawtooth in the frame budget.

use crate::layer::Layer;
use crate::time::Tick;

/// A system's declared contract with the engine.
#[derive(Clone, Debug)]
pub struct SystemDesc {
    pub name: &'static str,
    /// The layer this system is responsible for and writes to.
    pub writes: Layer,
    /// Every layer this system reads. Checked against `writes` at registration.
    pub reads: &'static [Layer],
    /// Run once every `period` ticks. Must be >= 1.
    pub period: u64,
    /// Offset within the period, to stagger equal-period systems. `< period`.
    pub phase: u64,
}

impl SystemDesc {
    pub fn new(name: &'static str, writes: Layer) -> Self {
        SystemDesc { name, writes, reads: &[], period: 1, phase: 0 }
    }

    pub fn reads(mut self, reads: &'static [Layer]) -> Self {
        self.reads = reads;
        self
    }

    pub fn every(mut self, period: u64) -> Self {
        self.period = period;
        self
    }

    pub fn phase(mut self, phase: u64) -> Self {
        self.phase = phase;
        self
    }
}

/// A unit of simulation.
///
/// Generic over the world type so that `au-core` stays ignorant of what a
/// world actually contains. The core knows about *time* and *ordering*; it
/// must never know about planets.
pub trait System<W> {
    fn run(&mut self, world: &mut W);
}

/// Blanket impl so a plain closure can be a system. Handy for tests and for
/// small systems that don't deserve a struct.
impl<W, F: FnMut(&mut W)> System<W> for F {
    fn run(&mut self, world: &mut W) {
        self(world)
    }
}

struct Entry<W> {
    desc: SystemDesc,
    order: usize,
    system: Box<dyn System<W>>,
    /// Diagnostics. Cheap, and the alternative is guessing.
    runs: u64,
}

/// Why a system was refused at registration.
#[derive(Debug, PartialEq, Eq)]
pub enum ScheduleError {
    /// A system tried to read a layer above the one it writes. This is the
    /// Golden Rule violation, caught before the universe ever starts.
    ReadsAbove { system: &'static str, writes: Layer, reads: Layer },
    ZeroPeriod { system: &'static str },
    PhaseOutOfRange { system: &'static str, period: u64, phase: u64 },
    DuplicateName { system: &'static str },
}

impl std::fmt::Display for ScheduleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ScheduleError::ReadsAbove { system, writes, reads } => write!(
                f,
                "LAYER VIOLATION: system '{}' writes {} but reads {}. \
                 A system may not read a layer above its own — that is the engine \
                 asking for content instead of simulating causes. \
                 Move the logic down, or move the system up.",
                system,
                writes.name(),
                reads.name()
            ),
            ScheduleError::ZeroPeriod { system } => {
                write!(f, "system '{}' has period 0; must be >= 1", system)
            }
            ScheduleError::PhaseOutOfRange { system, period, phase } => write!(
                f,
                "system '{}' has phase {} >= period {}; it would never run",
                system, phase, period
            ),
            ScheduleError::DuplicateName { system } => {
                write!(f, "system '{}' is registered twice", system)
            }
        }
    }
}

impl std::error::Error for ScheduleError {}

/// Dispatches systems by tick, in layer order.
pub struct Scheduler<W> {
    entries: Vec<Entry<W>>,
    sorted: bool,
}

impl<W> Default for Scheduler<W> {
    fn default() -> Self {
        Scheduler { entries: Vec::new(), sorted: true }
    }
}

impl<W> Scheduler<W> {
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a system. Rejects layer violations before anything runs.
    pub fn register(
        &mut self,
        desc: SystemDesc,
        system: impl System<W> + 'static,
    ) -> Result<(), ScheduleError> {
        if desc.period == 0 {
            return Err(ScheduleError::ZeroPeriod { system: desc.name });
        }
        if desc.phase >= desc.period {
            return Err(ScheduleError::PhaseOutOfRange {
                system: desc.name,
                period: desc.period,
                phase: desc.phase,
            });
        }
        if self.entries.iter().any(|e| e.desc.name == desc.name) {
            return Err(ScheduleError::DuplicateName { system: desc.name });
        }
        for &r in desc.reads {
            if !desc.writes.may_read(r) {
                return Err(ScheduleError::ReadsAbove {
                    system: desc.name,
                    writes: desc.writes,
                    reads: r,
                });
            }
        }
        let order = self.entries.len();
        self.entries.push(Entry { desc, order, system: Box::new(system), runs: 0 });
        self.sorted = false;
        Ok(())
    }

    /// Panicking form. Use in `main`/setup, where a layer violation should stop
    /// the world.
    pub fn must_register(&mut self, desc: SystemDesc, system: impl System<W> + 'static) {
        if let Err(e) = self.register(desc, system) {
            panic!("{}", e);
        }
    }

    fn ensure_sorted(&mut self) {
        if !self.sorted {
            // Bottom-up, then registration order within a layer.
            // Never by name, never by hash: both are silent determinism bugs.
            self.entries.sort_by_key(|e| (e.desc.writes, e.order));
            self.sorted = true;
        }
    }

    /// Run every system due on `tick`, in layer order.
    pub fn run_tick(&mut self, tick: Tick, world: &mut W) {
        self.ensure_sorted();
        for e in self.entries.iter_mut() {
            if tick.0 % e.desc.period == e.desc.phase {
                e.system.run(world);
                e.runs += 1;
            }
        }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// (name, layer, period, times run). For the debugger and for telemetry.
    pub fn report(&self) -> Vec<(&'static str, Layer, u64, u64)> {
        self.entries
            .iter()
            .map(|e| (e.desc.name, e.desc.writes, e.desc.period, e.runs))
            .collect()
    }
}
