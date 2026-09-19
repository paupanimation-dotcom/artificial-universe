//! # au-sim — where the layers are assembled
//!
//! This crate owns the `World` and the tick loop, and it is deliberately thin.
//! Layers are added by *registering systems*, not by editing the loop. If this
//! crate starts to grow logic, the architecture is leaking.

pub mod autocatalysis;
pub mod chem_columns;
pub mod chemistry;
pub mod config;
pub mod observe;
pub mod physics_columns;
pub mod protocell;
pub mod planet_gen;
pub mod sim;
pub mod systems;
pub mod world;

pub use config::{Config, ConfigError, DEFAULT_CONFIG};
pub use sim::Simulation;
pub use world::World;

/// A simulation from the built-in defaults, with an optional seed override.
pub fn boot(seed: Option<u64>) -> Simulation {
    let mut cfg = Config::parse(DEFAULT_CONFIG).expect("built-in config must parse");
    if let Some(s) = seed {
        cfg.set("world.seed", &s.to_string());
    }
    Simulation::new(cfg)
}
