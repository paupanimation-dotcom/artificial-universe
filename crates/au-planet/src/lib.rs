//! # au-planet — planets as physics, not painting
//!
//! A planet here is not a heightmap. It is a small set of declared physical
//! parameters — mass, radius, orbital distance, stellar luminosity, composition —
//! and everything at the largest scale *follows* from them by the physics already
//! built in Phase 2. Gravity from Newton, temperature from Stefan–Boltzmann,
//! ocean-or-ice-or-vapour from water's own phase diagram. Nobody paints an ocean;
//! the star's warmth and water's boiling point decide there is one.
//!
//! [`planet`] holds the parameter → state derivation, which is the whole of the
//! scientific claim: the habitable zone is not a band we drew, it is the orbital
//! range over which the phase model returns liquid water.
//!
//! [`terrain`] holds the one honestly-declared exception: relief detail, drawn
//! deterministically from the seed *within* the physical envelope — with the
//! licence and its limits documented where the dice are rolled.
//!
//! The chunk generator that assembles these into stored matter lives in `au-sim`
//! (as chemistry's did), because that is where planets meet columns; this crate
//! stays pure physics and pure functions, with no storage dependency at all.

pub mod planet;
pub mod terrain;

pub use planet::{
    Planet, PlanetParams, SurfaceState, AU, BIG_G, REFERENCE_PRESSURE, SOLAR_LUMINOSITY,
};
pub use terrain::Terrain;
