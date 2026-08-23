ARTIFICIAL UNIVERSE — PHASE 5a (OPEN-ENDED CHEMISTRY) — INCREMENTAL DELTA
═════════════════════════════════════════════════════════════════════════

PREREQUISITE: apply Phase 3, 3b, 4, and 4b deltas first, in order.
Overlay this archive onto your project folder, letting files overwrite.

WHAT'S IN HERE
──────────────
NEW:
    crates/au-chem/src/network.rs      the generative reaction network
    crates/au-chem/tests/network.rs    9 tests incl. the_engine_discovers_water
                                       and Boltzmann equilibrium on GENERATED
                                       reactions
    crates/au-sim/tests/genesis.rs     4 tests: discovery in the sim loop,
                                       determinism, resume-by-name, pool wall

CHANGED:
    crates/au-chem/src/lib.rs          exports
    crates/au-chem/src/molecule.rs     symmetry_classes() extracted from
                                       canonical() (canonical unchanged)
    crates/au-chem/src/reaction.rs     SpeciesRegistry::to_bytes/from_bytes
    crates/au-core/src/event.rs        CHEMISTRY_DISCOVERY event kind
    crates/au-sim/src/chemistry.rs     chem.open_ended config (seeds, pool,
                                       barrier, pre-exponential, max_atoms)
    crates/au-sim/src/world.rs         chem_registry on World: seeded at boot,
                                       folded into the world hash
    crates/au-sim/src/systems/chemistry.rs  the open mode: enumerate from what
                                       exists, intern, log, filter to the pool
    crates/au-sim/src/sim.rs           framed save: [len][snapshot][registry]
    crates/au-headless/src/main.rs     the `genesis` demo
    README.md, CHANGELOG.md

NOTE: the save format gained a frame; saves from earlier phases no longer
load. The au-data codec itself is untouched.

AFTER OVERLAYING
────────────────
    cargo test --release            # 128 tests (was 115)
    cargo run --release -- genesis  # H and O atoms; the engine discovers water

WHAT PHASE 5a PROVES
────────────────────
Reactions are DERIVED from molecular structure (four graph moves in exact
inverse pairs), molecules are DISCOVERED at runtime (interned, event-logged,
hashed, serialized), and the generated chemistry still lands on real
thermodynamics: Ea_f - Ea_r = dH holds by construction across the whole
closure, so K = exp(-dH/RT) — verified by running a reactor on purely
generated reactions and measuring Boltzmann equilibrium at two temperatures.

From H and O atoms alone the world discovers H2, OH, O2, H2O, HO2 — and the
demo surfaced kinetic trapping UNASKED: at 1500 K matter jams at hydroxyl
(58 molecules of water out of 2M possible); at 3000 K the barriers open and
99.86% of the oxygen anneals into water. Four moves, one activation identity,
and real chemistry's signature behaviour falls out.

HONEST LIMITS (stated in code and CHANGELOG)
────────────────────────────────────────────
One cell — species transport between cells is Phase 5b. Rings and concerted
substitutions excluded (documented pairs for later). Enthalpy-only equilibria
(no dS term yet). Pool = shelf space: discovery past it is recorded, never
stocked.
