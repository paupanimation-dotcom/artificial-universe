//! # au-core — the deterministic foundation
//!
//! Everything in this crate exists to defend one property: **the same seed
//! produces the same universe, on any machine, forever.**
//!
//! Not for its own sake. It is what makes the promises in PROJECT_VISION
//! physically possible:
//!
//! * "Rewind historical data" — you cannot rewind what you cannot reproduce.
//! * "Procedural reconstruction" (ARCHITECTURE §19) — a world too large to store
//!   must be regenerable from its seed, exactly.
//! * "Follow descendants, view family trees" — a lineage is a claim about the
//!   past, and it must still be true a billion ticks later.
//! * Debugging emergence at all. When no life appears after 400 million years,
//!   the only way to find out why is to replay it.
//!
//! This crate knows nothing about planets, molecules or organisms, and it never
//! will. It knows about time, order, identity and randomness.
//!
//! **Zero dependencies.** Deliberate. Every crate is a place where a version
//! bump can silently change an iteration order or a floating-point path.

pub mod event;
pub mod hash;
pub mod ids;
pub mod layer;
pub mod rng;
pub mod schedule;
pub mod time;

pub use event::{Event, EventLog};
pub use hash::{Hasher, WorldHash};
pub use ids::{EntityAllocator, EntityId};
pub use layer::Layer;
pub use rng::{Domain, Rng};
pub use schedule::{ScheduleError, Scheduler, System, SystemDesc};
pub use time::{SimClock, SimDuration, SimInstant, Tick};
