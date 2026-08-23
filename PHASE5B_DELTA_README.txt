Artificial Universe — Phase 5b delta (species transport: diffusion)
===================================================================

Apply AFTER the Phase 3, 3b, 4, 4b and 5a deltas, preserving paths
(overwrite when prompted):

    unrar x -o+ au-phase5b-delta.rar /path/to/project/

New files:
    crates/au-physics/src/diffuse.rs        Fick diffusion on integer counts:
                                            explicit + implicit (backward Euler,
                                            red-black GS), donor's-pocket
                                            limiter, sweeps_for_deep_time
    crates/au-physics/tests/diffuse.rs      10 tests: exact conservation at any
                                            convergence, Einstein σ²=2kt, walls,
                                            the deep-time n²-sweep regression
    crates/au-sim/tests/diffusion.rs        5 tests: Graham's law end-to-end,
                                            frozen host traps chemistry, deep
                                            time, discovered water spreading,
                                            save/resume mid-spread

Changed files:
    crates/au-physics/src/lib.rs            export the diffuse module
    crates/au-sim/src/chemistry.rs          chem.diffusion_m2_s (default 0=off)
    crates/au-sim/src/systems/chemistry.rs  Graham-scaled per-species rates from
                                            the discovered graphs; phase-derived
                                            mobility mask; diffusion pass between
                                            react and scatter; deep-time path on
                                            the SAME physics.implicit_conduction
                                            flag as heat
    crates/au-headless/src/main.rs          `au vent` — the two-act demo
    README.md, CHANGELOG.md                 status, 143 tests, the full story
                                            (both bugs documented honestly)

Verify:
    cargo test --release            # 143 tests
    cargo run --release -- vent     # water invented at one cell, found at the ends
