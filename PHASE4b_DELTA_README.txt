ARTIFICIAL UNIVERSE — PHASE 4b (IMPLICIT CONDUCTION) — INCREMENTAL DELTA
════════════════════════════════════════════════════════════════════════

PREREQUISITE: apply Phase 3, 3b, and 4 deltas first, in order.
Overlay this archive onto your project folder, letting files overwrite.

WHAT'S IN HERE
──────────────
NEW:
    crates/au-physics/src/implicit.rs     backward-Euler conduction solver
    crates/au-physics/tests/implicit.rs   12 tests: stability, agreement,
                                          conservation-under-bad-convergence

CHANGED:
    crates/au-physics/src/lib.rs          module + exports
    crates/au-physics/src/grid.rs         Boundary: implicit_conduction,
                                          conduction_iters, conduction_omega
    crates/au-physics/src/transport.rs    scratch, stable_dt exemption, routing,
                                          and the infinite-dt sentinel fix
    crates/au-sim/src/systems/physics.rs  config: physics.implicit_conduction etc.
    crates/au-headless/src/main.rs        the `deeptime` demo
    README.md, CHANGELOG.md

AFTER OVERLAYING
────────────────
    cargo test --release             # 115 tests (was 103)
    cargo run --release -- deeptime  # the dx² wall, and the chain coming off

WHAT PHASE 4b PROVES
────────────────────
The project's oldest limit is gone. An explicit diffusion scheme is stable only
while dt <= dx^2/(2*ndim*alpha) — halve the cell size and the timestep quarters,
so a million years of rock at metre resolution costs ~3e8 substeps for ONE
chunk. Backward Euler removes that bound entirely: measured 5000x fewer
substeps and ~800x less wall clock for the same span, with the advantage
growing without bound as the grid is refined.

Crucially, exactness survived. An iterative solver converges to a TOLERANCE,
but the engine's energy is an exact i128. The resolution: the SOLVER decides
how much energy should move (approximate); the TRANSPORT moves it as one
integer per face, subtracted from one side and added to the other (exact).
Conservation therefore does not depend on the solve converging — a test runs
the solver crippled at ONE sweep and the books still balance to exactly 0 uJ.

TWO REAL BUGS, BOTH CAUGHT BY TESTS (documented in code + CHANGELOG)
────────────────────────────────────────────────────────────────────
1. stable_dt returns infinity when nothing bounds the step. That used to mean
   only "inert grid", so an early return was safe. Implicit conduction gave it
   a second meaning and the shortcut silently skipped the whole physics step.
   Notably the CONSERVATION tests passed throughout — doing nothing conserves
   beautifully — and only the physics tests caught it.
2. Passing the pressure solver's omega=1.85 to the heat solve produced 4244 K
   in a box that started between 400 K and 1200 K: conserved, but physically
   impossible. The conduction update is a convex combination (positive weights
   summing to one), which cannot leave its starting range — unless you
   over-relax, which breaks exactly that property. conduction_omega now
   defaults to 1.0.

THE LIMIT THAT REMAINS
──────────────────────
Stability is unconditional; ACCURACY is not — a big step smears fast transients.
Radiation (T^4) and advection (dx/u) are still explicit and still bound the
step, but neither scales as 1/dx^2, so neither worsens under refinement.
Viscosity has the same dx^2 problem and is the next candidate; implicit.rs is
written to be reused for it. The cliff is gone; the hills remain.
