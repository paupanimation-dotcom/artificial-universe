//! Simulation systems, one module per layer as they are built.
//!
//! * **Universe** — `time_scale`: how fast the clock *wants* to run.
//! * **Physics** — `physics`: thermodynamics, and how fast the clock is
//!   *allowed* to run. Those two are in tension, deliberately, and physics wins.
//!
//! Layers 2..13 are empty. Not stubbed — empty. A stub that returns plausible
//! numbers is worse than nothing, because everything above it gets built against
//! a fiction and nobody finds out for a year.

pub mod chemistry;
pub mod physics;
pub mod time_scale;
pub mod protocell;
