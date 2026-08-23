ARTIFICIAL UNIVERSE — PHASE 4 (PLANET GENERATION) — INCREMENTAL DELTA
═════════════════════════════════════════════════════════════════════

PREREQUISITE: apply the Phase 3 and Phase 3b deltas first, in order.
Overlay this archive onto your project folder, letting files overwrite.

WHAT'S IN HERE
──────────────
NEW — the planet crate and its generator:
    crates/au-planet/Cargo.toml
    crates/au-planet/src/lib.rs         planets as physics, not painting
    crates/au-planet/src/planet.rs      parameters → gravity, temperature, ocean/ice/vapour
    crates/au-planet/src/terrain.rs     seeded relief — the licence, documented
    crates/au-planet/tests/validate.rs  10 tests incl. the emergent habitable zone
    crates/au-sim/src/planet_gen.rs     the chunk generator: seed decides where, physics decides what
    crates/au-sim/tests/planet.rs       5 integration tests incl. warm→ocean / cold→ice

CHANGED:
    Cargo.toml                          au-planet in the workspace
    crates/au-core/src/rng.rs           PLANET RNG domain (position-keyed relief)
    crates/au-sim/Cargo.toml, lib.rs    dependency + module
    crates/au-sim/src/sim.rs            planet generator auto-installed on boot AND resume
    crates/au-headless/*                the `planet` demo command
    README.md, CHANGELOG.md

AFTER OVERLAYING
────────────────
    cargo test --release           # 103 tests (was 84)
    cargo run --release -- planet  # Earth's numbers derived; the habitable zone;
                                   # the same seed as sea at 0.95 AU, ice at 1.6 AU

WHAT PHASE 4 PROVES
───────────────────
A world can declare a planet (mass, radius, orbit, star, albedo) and its chunks
fill themselves — the first matter in the project not placed by a fixture. The
gross state is DERIVED by the physics already built: gravity (Newton),
temperature (Stefan–Boltzmann — 255.5 K for Earth, the textbook number), and
ocean-vs-ice-vs-vapour by water's own phase model. The habitable zone emerges
as a contiguous band nobody drew, and it slides outward with brighter stars.
The headline test: the SAME generator at 0.95 AU makes an ocean and at 1.6 AU
makes an ice sheet — physics decided, never the generator. Determinism and
save/resume hold; resume now reinstalls the planet automatically.

Two honest limits the physics surfaced, both documented in code: a bare-rock
Earth freezes (the missing 33 K IS the greenhouse, not yet modelled), and Venus
would hold a hot liquid ocean without its runaway atmosphere (my first test
assumed otherwise; the engine failed it correctly).

DEFERRED (written at the code sites they constrain)
───────────────────────────────────────────────────
Circulation, weather, hydrology (need the implicit solver); greenhouse
radiative transfer; tectonics (relief becomes an OUTPUT once it exists);
partial-cell coastlines; seasons; derived surface pressure.
