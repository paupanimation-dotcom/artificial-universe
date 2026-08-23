ARTIFICIAL UNIVERSE — PHASE 3 (CHEMISTRY) — INCREMENTAL DELTA
════════════════════════════════════════════════════════════

This archive contains ONLY the files that are new or changed since the
Phase 2b delta. Overlay it onto your existing project folder, letting the
files overwrite in place. The folder structure is preserved.

WHAT'S IN HERE
──────────────
NEW — the whole chemistry crate:
    crates/au-chem/Cargo.toml
    crates/au-chem/src/lib.rs          crate overview & the honest framing
    crates/au-chem/src/element.rs      elements as exactly-conserved i128 counts
    crates/au-chem/src/molecule.rs     molecules as graphs + canonical form
    crates/au-chem/src/reaction.rs     species registry, reactions, bond-energy model
    crates/au-chem/src/reactor.rs      the reactor: chemistry coupled to Phase 2 heat
    crates/au-chem/tests/foundation.rs 9 tests: conservation + graph identity
    crates/au-chem/tests/reactions.rs  9 tests: kinetics, equilibrium, van 't Hoff

CHANGED:
    Cargo.toml                         added au-chem to workspace members
    crates/au-headless/Cargo.toml      added au-chem dependency
    crates/au-headless/src/main.rs     added the `react` demo command
    crates/au-physics/src/fluid.rs     removed a no-op drop() (warning cleanup only)
    README.md                          status → Phase 3, added `react` to commands
    CHANGELOG.md                       prepended the Phase 3 section

AFTER OVERLAYING
────────────────
    cargo test --release          # 79 tests (was 61)
    cargo run --release -- react  # chemistry heats its own cell; equilibrium shifts

WHAT PHASE 3 PROVES
───────────────────
The engine reproduces van 't Hoff — an exothermic reaction's equilibrium
constant falls with temperature, emerging from nothing but two Arrhenius
rates. Exothermic reactions heat their own cell through the SAME energy field
Phase 2 conducts and radiates, with conservation exact to the microjoule.
"No hardcoded molecules" is real: water is recognised as a graph, never a
table entry.

DEFERRED TO PHASE 3b
────────────────────
Wiring chemistry into au-sim's scheduler (per-cell species storage, a
chemistry system on the active set, snapshot-format extension). The engine is
validated in isolation and via the demo; integration is deferred so it can be
done carefully. See the CHANGELOG for details.
