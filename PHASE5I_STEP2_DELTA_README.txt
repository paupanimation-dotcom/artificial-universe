Artificial Universe — Phase 5i step 2 delta
===========================================

PROTOCELLS: COMPOSITION DETERMINES REPLICATION RATE.

9 files. 4 new, 5 modified. Extract over your project and overwrite.
Apply Phase 5h and 5i step 1 first.

  NEW
    crates/au-sim/src/protocell.rs          Protocell, ProtocellStore, lineage
    crates/au-sim/src/systems/protocell.rs  nucleation, uptake, chemistry, fission
    crates/au-sim/tests/protocell.rs        7 tests

  MODIFIED
    CHANGELOG.md                            Phase 5i step 2 entry
    README.md                               status + Phase 5i step 2 section
    crates/au-sim/src/lib.rs                module registration
    crates/au-sim/src/systems/mod.rs        system registration
    crates/au-sim/src/world.rs              World::protocells - hashed, restored
    crates/au-sim/src/sim.rs                life trailer, system wiring

VERIFIED
  cargo test --release  ->  237 passed, 0 failed  (was 230), 37 binaries
  Built from a clean tree with rustc 1.75; zero third-party dependencies.

THE RESULT
  A bag imports substrate, converts it to surfactant inside itself, and the
  surfactant cannot leave because it IS the membrane. Area ratchets up, the
  reduced volume crosses 1/sqrt(2), the bag divides, and both daughters
  inherit the composition that caused it. A world whose chemistry cannot
  build membrane never divides. Nothing grants membranes an advantage.

OFF BY DEFAULT
  life.protocells = false. A world that does not ask for them is unchanged.
  Optional keys: life.bilayer_thickness_m, life.lysis_strain,
  life.min_radius_thicknesses.

SAVE FORMAT
  Frames gain a second optional trailer ("AULIFE01") carrying protocells and
  the entity allocator, written after the atom trailer and peeled before it.
  Pre-5i saves load with no protocells. No existing format was revised.

LOGGED, NOT SOLVED
  - Bond-derived reaction enthalpy for C5H12 + H2O -> C5H12O + H2 comes out
    ~11 orders of magnitude above the bond balance. Phase 3b question; the
    test declares enthalpy_j rather than silently measuring it.
  - clock.rescale.enabled = false did not stop rescaling (dt = 3.6e-4).
  - Declared reactions only inside bags; open-ended generation needs the
    reaction cache shared between two systems.
