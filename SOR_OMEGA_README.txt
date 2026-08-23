Artificial Universe — the relaxation factor is now derived, not declared
=======================================================================

4 files. Extract over your project and overwrite.
Apply the implicit-viscosity delta first (this builds on it).

VERIFIED
  au-physics: 51 passed, 0 failed, 1 ignored (was 49). Ra_c still found by
  bisection. NOT re-run across the whole workspace.

THE DEFECT
  Optimal SOR for a Laplace problem is 2/(1+sin(pi/n)) - a property of the
  GRID, not a constant. sor_omega was declared 1.85, which is optimal for a
  41-cell grid: the tank Phase 2b validated Ra_c on. Wrong for every grid
  since.

  Measured, 12-cell grid, 20 sweeps:
    omega 1.00  ->  83.4% divergence
    omega 1.40  ->  15.4%
    omega 1.59  ->   1.8%   <- derived
    omega 1.85  -> 185.1%   <- the declared constant
    omega 1.95  -> 168.3%

  The declared value was WORSE THAN NO OVER-RELAXATION AT ALL. Past the
  optimum SOR oscillates rather than converging slowly.

  sor_omega now defaults to 0 = derive. Nonzero is still honoured.

ALSO MEASURED: THE PROJECTION IS THE COST
  20 ticks, 24x12 cell: 4.09 s at 60 sweeps, 0.20 s at 15. But 15 sweeps
  gives 200% divergence - a compressible fluid, whose drifted density
  produces buoyancy indistinguishable from real convection. 60 sweeps is
  near the edge of honesty, NOT padding. Do not cut it.

WHAT THIS DID NOT BUY
  The convection cell got slower (5.34s vs 4.09s for 20 ticks). Hypothesis:
  an honest projection lets the flow actually convect, so velocities rise
  and the Courant bound tightens. UNTESTED. Correctness improved; measured
  throughput did not.

WHERE SPEED HAS TO COME FROM
  Multigrid. Single-grid Gauss-Seidel converges at a rate set by the grid
  (~cos^2(pi/n) per sweep) - the same finding Phase 5b recorded for deep-time
  diffusion. One multigrid would serve conduction, species transport and the
  projection together.
