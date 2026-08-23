Artificial Universe — complete project, assembled and verified
==============================================================

This is THE tree. Not a delta, not an overlay. Delete whatever you had and
use this. Nothing needs to be applied on top of anything.

VERIFIED, on this exact tree, just now:
    cargo build --release      clean
    cargo test  --release      252 passed, 0 failed, 6 ignored, 38 binaries

AUDITED
  Every feature built in this project was checked by name against this tree
  after assembly. Seven were missing (the lysis/washout counter split, the
  turnover tests, the two diagnostic probes, and the cascade analysis) - all
  from work done after the last delta was packaged, lost to a container
  reset. They have been reconstructed and the suite re-run: 250 -> 252.

HOW IT WAS BUILT
  Base: au-full-current.tar.gz (complete tree through Phase 5g, 71 source
  files) with these applied in order:
      5h open system -> 5i step 1 vesicles -> 5i step 2 protocells ->
      enthalpy units fix -> observation -> protocell drain -> uptake
      limiter -> bounded nucleation -> implicit viscosity -> derived
      sor_omega -> Brownian motion
  81 source files. The copy you sent was missing au-chem and au-planet
  entirely and stopped at the drain.

BUILD
    apt-get install -y rustc cargo        (1.75; zero third-party deps)
    cargo build --release
    cargo test --release
    ./target/release/au watch --ticks 3000 --out world.html

  Other demos: physics convect ra react react-sim planet deeptime genesis
  vent raf

WHAT IS IN THE CONTEXT PANEL, NOT HERE
  PROJECT_VISION.md, ARCHITECTURE.md, DEVELOPMENT_ROADMAP.md are in this
  tree for reference, but the *_for_project_context.md variants are the ones
  to upload to the Claude project's Context panel. Newest wins.

WHERE THE SCIENCE ACTUALLY STANDS
  Phases 1-5i are built and validated against anchors that predate the code:
  Fourier, Stefan, Clausius-Clapeyron, Rayleigh's Ra_c = 657.5 found by
  bisection, Graham's law, Hordijk-Steel RAF, the CMC, Fick's steady
  profile, the two-sphere bound v = 1/sqrt(2), binomial 1/(2*sqrt(n)),
  Stokes-Einstein.

  Phase 5j is in progress: protocells nucleate, feed, run their own
  chemistry, divide, record lineage, wander by Brownian motion, and wash out
  through an open port.

THE OPEN PROBLEM - READ THIS BEFORE BUILDING ON PROTOCELLS
  Nucleation is correct (traced: captures osmotic balance, ratio 1.000).
  But each nucleation removes ~1675 osmolytes from its grid cell, so c_ext
  falls, so every bag already floating there is above the new external
  concentration and sheds solute to match. Its osmotic volume shrinks at
  fixed membrane area, v crosses 1/sqrt(2), and it fissions.

  NUCLEATION TRIGGERS A FISSION WAVE ACROSS THE STANDING POPULATION, and
  each wave doubles it. Measured: 46,497 divisions from 1,020 nucleations in
  a world with NO REACTIONS DECLARED.

  Consequence: DIVISION COUNTS ARE NOT EVIDENCE OF GROWTH. The Phase 5i
  step 2 claim that composition determines replication rate must be
  re-earned against this, not assumed to have survived.

  Part is real physics - draw solute from a medium and vesicles shrink. What
  is not real is the coupling: bags are made from the same finite cell they
  float in, so their own creation squeezes their neighbours. A real medium
  is not measurably depleted by one vesicle forming in it. Candidate fixes:
  a reservoir holding solute fixed (as Phase 5h holds substrate), or
  nucleation drawing osmolytes from outside the cell budget. NOT YET CHOSEN.

ALSO NOT DONE, HONESTLY
  - 6 ignored tests, each with its reason at the test: two convection ones
    are too slow (the pressure projection dominates), one measures nothing
    on a grid where viscosity does not bind, one is a diagnostic probe.
  - Multigrid is the outstanding physics work. Single-grid Gauss-Seidel
    converges at a rate set by the grid; one multigrid would serve
    conduction, species transport and the pressure projection together.
