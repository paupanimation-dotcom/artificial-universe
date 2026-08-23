Artificial Universe — Phase 5h delta
====================================

THE OPEN SYSTEM: matter can now cross a boundary.

10 files. 2 new, 8 modified. Extract over your project folder and overwrite.
The other 40-odd files in the tree are untouched.

  NEW
    crates/au-chem/src/reservoir.rs        AtomLedger, Port, port_delta
    crates/au-sim/tests/openworld.rs       11 integration tests

  MODIFIED
    CHANGELOG.md                           5h entry + 5e/5f/5g reconstructed
    README.md                              status, stack, Phase 5h section
    crates/au-chem/src/lib.rs              module registration + exports
    crates/au-sim/src/chemistry.rs         port config parsing and validation
    crates/au-sim/src/sim.rs               atom-ledger trailer in the save frame
    crates/au-sim/src/world.rs             World::atoms — hashed, restored
    crates/au-sim/src/systems/chemistry.rs ports applied per tick, gross booking
    crates/au-sim/tests/membrane.rs        +2 tests: a membrane in a flow

VERIFIED
  cargo test --release  ->  215 passed, 0 failed  (was 196)
  Built from a clean tree with rustc 1.75; zero third-party dependencies.

A NOTE ON THE SAVE FRAME
  Saves gain an optional trailer at the very end, identified by the magic
  "AUATOMS1", carrying the per-element atom ledger. A save written before this
  phase has no trailer and loads with an empty ledger, which is the correct
  reading of a world in which nothing could cross a boundary. No existing
  format was revised.

DEFAULTS ARE UNCHANGED
  A world that declares no ports runs byte-identically to Phase 5g. The two
  new config keys are:
      chem.reservoir.source = "0:1@200000, 2@200000; 5:1@50"
      chem.reservoir.sink   = "40; 41"
