Artificial Universe — reaction enthalpy units fix
=================================================

8 files. Extract over your project and overwrite.
Also included, OUTSIDE the au/ folder:
  ARCHITECTURE_for_project_context.md
    -> upload this to the project's Context panel, replacing ARCHITECTURE.md.
       It gains a section "Units, resolution, and things that cannot be
       represented" recording every trap that cost a session.

VERIFIED
  cargo test --release  ->  237 passed, 0 failed, 37 binaries.

THE DEFECT
  base_per_order is J/mol, so reaction_enthalpy was molar — but its doc said
  "per event" and two of three callers believed it. Declared reactions without
  an explicit enthalpy_j delivered 6.022e23x too much heat.

THE FIX
  reaction_enthalpy_molar() returns f64 J/mol (not a quantised Energy), so the
  conversion cannot be skipped. Both wrong callers now divide by AVOGADRO.
  network.rs was correct all along and is only renamed.

DO NOT "SIMPLIFY" THIS BY DIVIDING INSIDE THE FUNCTION
  That was tried. It breaks generated_reactions_reach_boltzmann_equilibrium
  (measured K -> 0 vs expected 8.9554e-29) because network.rs then divides
  twice. The anchor test is right.

TRAPS THIS UNCOVERED
  - Reaction::enthalpy (Energy) is microjoule-quantised. A correct per-event
    enthalpy (~5e-19 J) rounds to ZERO there. Read enthalpy_j.
  - Below ~2e12 reaction events the thermal field cannot resolve chemistry.
    Any "chemistry heats its cell" test needs populations above that.
  - activation_j is MOLAR (meets R in Arrhenius); enthalpy_j is PER EVENT
    (meets the thermal field). Do not mix them.
