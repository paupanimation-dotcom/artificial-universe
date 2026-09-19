//! # au-physics — thermodynamics and transport
//!
//! ## What physics means here
//!
//! Not rigid bodies. Not collisions. Not ragdolls.
//!
//! Look at what the layers above this one actually need. Chemistry needs
//! temperature, pressure and available energy. Planets need heat flow. Climate
//! *is* heat transport. And the origin of life needs one thing above all others:
//! a **persistent free-energy gradient**, because life is a dissipative structure
//! — it exists because the star is hot and the sky is cold, because the vent is
//! hot and the ocean is cold. It is what grows in the gap between them.
//!
//! Nobody in this stack needs a falling crate. Everybody needs *how much energy
//! is here, and where is it going*.
//!
//! So Phase 2 is thermodynamics, and building a rigid-body engine instead would
//! have been building the wrong thing beautifully.
//!
//! ## The one thing that must not be wrong
//!
//! Conservation. **Evolution is an adversarial optimiser whose search space
//! includes the bugs in your physics.** Leave any path by which energy can be
//! created and something will eventually evolve to walk it — not because it is
//! clever, but because natural selection is an exhaustive search and free energy
//! is the most rewarding thing in any possible fitness landscape.
//!
//! Hence: energy and mass are **exact integers**, transport is **symmetric
//! flux**, and conservation is an identity rather than an aspiration. See
//! [`quantity`].
//!
//! ## What this crate does not know about
//!
//! Worlds, chunks, planets, organisms. It operates on a [`grid::GridView`] —
//! three slices and a shape — and has no idea what it is a grid *of*. `au-sim`
//! is what connects it to a universe. That separation is not tidiness; it is what
//! lets the physics be tested against analytic solutions with no simulation
//! running at all.

pub mod fluid;
pub mod diffuse;
pub mod grid;
pub mod implicit;
pub mod material;
pub mod quantity;
pub mod state;
pub mod transport;

pub use fluid::{FluidScratch, FluidState};
pub use grid::{Boundary, FluidView, GridSpec, GridView, Wall};
pub use material::{Material, MaterialId, MaterialRegistry, Phase, P_REF, SIGMA};
pub use quantity::{Energy, Ledger, Mass, Momentum, P_PER_KG_M_S, AG_PER_KG, PJ_PER_J};
pub use state::{derive, energy_for_temperature, CellState};
pub use diffuse::{
    diffuse_explicit, diffuse_explicit_perm, diffuse_implicit, diffuse_implicit_perm,
    sweeps_for_deep_time, DiffuseScratch,
};
pub use implicit::{conduct_implicit, ImplicitScratch};
pub use transport::{step, Scratch, StepReport};
