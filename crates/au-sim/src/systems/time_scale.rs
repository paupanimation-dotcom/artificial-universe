//! Deep-time acceleration.
//!
//! ARCHITECTURE §3 asks for "long-term simulation acceleration", and it is not
//! optional: chemistry is interesting at the scale of seconds, evolution at the
//! scale of millions of years. The gap is ~10^14. No machine will ever run 10^14
//! ticks. So the clock coarsens as the universe ages — seconds per tick doubles
//! on a schedule until it reaches a ceiling.
//!
//! # This is a lossy trade, and it should be argued about
//!
//! Coarsening time does not merely make the universe run *faster*. It makes it
//! run *differently*. Fast events — a lightning strike, a reaction that
//! completes in milliseconds — stop being resolvable once a tick is a year
//! long. A world that fast-forwards through its Hadean is a different world,
//! not the same world seen sooner.
//!
//! The honest long-term answer is not one global clock but **per-region time**:
//! a tide pool where prebiotic chemistry is doing something interesting keeps a
//! fine clock, while the empty ocean floor runs coarse. The [`au_core::Scheduler`]
//! is already built for that — a system registers a *period*, and different
//! regions can register the same layer at different periods.
//!
//! We are not building that today. A single global scale is the right Phase 1
//! answer, and it is deliberately data-driven so that changing this policy later
//! costs a config edit, not an architecture.

use au_core::event::kind;
use au_core::layer::Layer;
use au_core::schedule::{System, SystemDesc};
use au_core::time::SimDuration;

use crate::world::World;

pub struct TimeScalePolicy {
    every_ticks: u64,
    factor: u64,
    max_secs: u64,
}

impl TimeScalePolicy {
    pub fn from_config(world: &World) -> Option<(SystemDesc, TimeScalePolicy)> {
        if !world.config.bool_or("clock.rescale.enabled", false) {
            return None;
        }
        let every_ticks = world.config.u64_or("clock.rescale.every_ticks", 65_536).max(1);
        let factor = world.config.u64_or("clock.rescale.factor", 2).max(1);
        let max_secs = world.config.u64_or("clock.rescale.max_seconds_per_tick", 31_557_600);

        let desc = SystemDesc::new("universe.time_scale", Layer::Universe)
            .reads(&[Layer::Universe])
            .every(every_ticks)
            // Phase 0 within the period would fire on tick 0, before the
            // universe has aged at all. Fire at the *end* of the first period.
            .phase(every_ticks - 1);

        Some((desc, TimeScalePolicy { every_ticks, factor, max_secs }))
    }

    pub fn every_ticks(&self) -> u64 {
        self.every_ticks
    }
}

impl System<World> for TimeScalePolicy {
    fn run(&mut self, world: &mut World) {
        let current = world.clock.scale();
        if current.secs >= self.max_secs {
            return;
        }
        // Integer multiply: the clock must stay bit-exact. An f64 factor here
        // would introduce the one thing time.rs exists to prevent.
        let mut next = current.mul_u64(self.factor);
        if next.secs >= self.max_secs {
            next = SimDuration::from_secs(self.max_secs);
        }
        if next == current {
            return;
        }
        world.clock.set_scale(next);
        // The change is part of the causal record. When someone asks in Phase 9
        // why a lineage vanished in a single tick, this is the first thing to
        // check.
        world.emit(
            Layer::Universe,
            kind::TIME_SCALE_CHANGED,
            next.secs,
            next.frac,
            0,
        );
    }
}
