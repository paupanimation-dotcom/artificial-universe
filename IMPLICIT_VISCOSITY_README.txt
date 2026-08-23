Artificial Universe — implicit viscosity
========================================

7 files. Extract over your project and overwrite.
OFF BY DEFAULT: physics.viscous_iters = 0 reproduces every earlier phase.

VERIFIED
  au-physics: 49 passed, 0 failed, 1 ignored (was 47).
  NOT re-run across the whole workspace after the last edit - the au-sim
  protocell tests were green before it and the change there is a config line
  plus comments, but the full number was not observed.

WHAT IT IS
  Backward Euler on the viscous term, Gauss-Seidel per velocity component.
  The last member of the trio Phase 4b began.

  Integration is one word: viscous_solve() proposes a velocity field and
  moves nothing; the EXISTING symmetric-flux loop then reads that field
  instead of the current one. Solver proposes, transport disposes - third
  time this project has used that split.

  omega = 1.0, never over-relaxed. The update is a convex combination and
  over-relaxation forfeits the maximum principle.

THE LOAD-BEARING TEST
  A solve crippled to ONE sweep, run 500 steps: momentum balances exactly
  against the impulse ledger. Conservation does not depend on convergence,
  because both halves of a symmetric integer transfer are the same integer.

  Do NOT rewrite this as dp = m*(u_new - u_old) per cell. That is
  conservative only on a perfect solve.

THE HONEST MEASUREMENT
  Built to unblock the protocell convection cell. Measured there:

    20 ticks, 24x12 cell:  explicit 7.8s -> implicit 4.6s   (1.7x)

  1.7x, NOT orders of magnitude. Viscosity was never the dominant cost -
  an assumption the roadmap and I both made without measuring. Protocell
  mobility is still unverified for the same performance reason.

  What remains: the advective Courant bound (dt <= dx/u, which no implicit
  VISCOUS scheme removes - advection is hyperbolic) and the pressure
  projection's 60 sweeps per substep. Measure both before assuming
  fine-grid flow is affordable.

ONE TEST IGNORED, AND WHY
  implicit_viscosity_takes_the_viscous_term_off_the_clock: THE TEST IS
  WRONG, NOT THE CODE. The Rayleigh-Benard tank is calibrated where a tick
  already fits the stable bound, so both paths take 1 substep and the
  comparison measures nothing. Needs a case where viscosity actually binds.
