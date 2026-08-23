ARTIFICIAL UNIVERSE — PHASE 3b (CHEMISTRY IN THE SIM LOOP) — INCREMENTAL DELTA
══════════════════════════════════════════════════════════════════════════════

PREREQUISITE: apply the Phase 3 delta FIRST. This delta builds on the au-chem
crate that Phase 3 introduced. Overlay Phase 3, then overlay this on top.

This archive contains ONLY the files new or changed since the Phase 3 delta.
Overlay onto your project folder, letting files overwrite in place.

WHAT'S IN HERE
──────────────
NEW:
    crates/au-sim/src/chem_columns.rs       per-cell species columns (0x2xxx range)
    crates/au-sim/src/chemistry.rs          the world's declared chemistry, from config
    crates/au-sim/src/systems/chemistry.rs  the chemistry system: reactor in the loop
    crates/au-sim/tests/chemistry.rs        5 integration tests

CHANGED:
    crates/au-core/src/event.rs             added CHEMISTRY_HEAT event kind
    crates/au-sim/Cargo.toml                depends on au-chem
    crates/au-sim/src/lib.rs                new modules
    crates/au-sim/src/world.rs              register chemistry columns at boot
    crates/au-sim/src/sim.rs                register chemistry system (new + load paths)
    crates/au-sim/src/systems/mod.rs        expose the chemistry module
    crates/au-headless/src/main.rs          the `react-sim` demo command
    README.md, CHANGELOG.md

AFTER OVERLAYING
────────────────
    cargo test --release             # 84 tests (was 79)
    cargo run --release -- react-sim # chemistry in the sim loop; save/resume identical

WHAT PHASE 3b PROVES
────────────────────
Chemistry now runs INSIDE the simulation, in the scheduler, on the active set,
writing the SAME energy field physics conducts and radiates. An exothermic
reaction heats the cell physics then cools; the temperature feeds back into the
reaction rates. Five integration tests prove: heat enters the shared field,
atoms are conserved, the world stays deterministic (same seed → same hash), and
save/resume is bit-identical with chemistry active.

THE SCOPE (honest)
──────────────────
A world DECLARES its chemistry vocabulary (elements, molecules, reactions) in
config at boot — like the material table. This fixes what molecules CAN exist,
never the outcome (which reactions fire, how hot, whether it ignites — all
emergent). Open-ended molecular discovery is deferred to Phase 5+, where the
species registry becomes serialized state. The layers above won't change.
