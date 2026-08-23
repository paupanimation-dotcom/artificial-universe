Artificial Universe — material-bounded nucleation, and a named error class
==========================================================================

3 files + 1 doc. Extract au/ over your project and overwrite.
ARCHITECTURE_for_project_context.md goes in the project's Context panel,
replacing ARCHITECTURE.md.

VERIFIED
  Last full clean run: 246 passed, 0 failed, 1 ignored, 38 binaries.
  The nucleation change and its effect (1 -> 64 divisions) were measured
  directly via `au watch`.

  NOT re-confirmed after the final #[ignore] markers were applied - the
  container reset twice while running the suite. The markers are additive
  and cannot change a pass count, but they were not observed green.

WHAT CHANGED
  Nucleation closed ONE bag per cell per tick - the shape of a for loop, not
  a physical rate. The medium held 42,000,000 assembled surfactant against a
  bag's 12,567 (material for 3,300 vesicles) and dribbled it out over 4,200
  ticks, so bags appeared long after their substrate was gone.

  Now every bag the surface supports closes, capped at 64 per cell per tick.
  That cap is DECLARED and bounds work, not physics: leftover surface closes
  next tick.

    before: 4130 nucleated,  1 divided, 0 lysed
    after:  4130 nucleated, 64 divided, 0 lysed

THE NEW CEILING IS COMPETITION
  Exactly 64 divisions across a tenfold range of starting surfactant: the
  first tick's bags divide, then 4,000 bags share one cell's substrate and
  the population freezes. That makes WASHOUT the binding constraint, which
  needs mobility, which needs implicit viscosity.

MOBILITY: UNVERIFIED AFTER FOUR ATTEMPTS
  Each attempt found something real; none produced a moving bag.
   1. a uniform current cannot persist in a sealed tube (incompressibility)
   2. physics.fluid.enabled defaults false - the column existing is not the
      solver running
   3. the physics system gathers-solves-scatters, so an initial condition is
      overwritten. A FLOW CANNOT BE DECLARED, IT MUST BE DRIVEN (buoyancy or
      a pressure boundary), the way au convect has since Phase 2b
   4. a driven convection cell works but is far too slow to run in the suite

  IMPLICIT VISCOSITY IS THE BLOCKER and is the next thing to build.

THE DOC
  ARCHITECTURE.md gains "A number means nothing without its container" - one
  error class this project has hit SEVEN times since Phase 3 and documented
  as seven separate traps. Instances, never the pattern, which is why each
  looked new.
