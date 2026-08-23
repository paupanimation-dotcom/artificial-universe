Artificial Universe — Phase 5c delta (catalysis + autocatalytic sets)
=====================================================================

Apply AFTER the Phase 3, 3b, 4, 4b, 5a and 5b deltas, preserving paths
(overwrite when prompted):

    unrar x -o+ au-phase5c-delta.rar /path/to/project/

New files:
    crates/au-chem/src/catalysis.rs     derived catalysis: the template model,
                                        absolute barrier reduction that leaves
                                        every equilibrium exactly untouched
    crates/au-chem/src/raf.rs           RAF detection (Hordijk & Steel 2004):
                                        closure, maximal_raf, self_producing
    crates/au-chem/tests/raf.rs         14 tests: the structural impossibility
                                        proof, the thermodynamic law, the
                                        algorithm against known answers, the
                                        phase transition
    crates/au-sim/src/autocatalysis.rs  survey(): a read-only lens on a world's
                                        discovered network (enumerates on a
                                        CLONE — looking cannot change the world)
    crates/au-sim/tests/autocatalysis.rs 7 tests: opt-in, conservation under
                                        catalysis, acceleration, the survey's
                                        non-interference, resume

Changed files:
    crates/au-chem/src/lib.rs           export catalysis + raf
    crates/au-sim/src/lib.rs            register the autocatalysis module
    crates/au-sim/src/chemistry.rs      chem.catalysis (default OFF) and its
                                        dials; refuses declared reactions
    crates/au-sim/src/systems/chemistry.rs  catalysed variants appended to the
                                        reaction set, capped deterministically
    crates/au-headless/src/main.rs      `au raf` — the three-act demo
    README.md, CHANGELOG.md             status, 164 tests, the full story
                                        (all three errors documented)

Verify:
    cargo test --release          # 164 tests
    cargo run --release -- raf    # C + O + CO -> 2 CO, and the collapse
