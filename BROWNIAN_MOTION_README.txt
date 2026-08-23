Artificial Universe — protocells move (Brownian motion)
=======================================================

3 files. Extract over your project and overwrite.

VERIFIED
  cargo test --release -> 250 passed, 0 failed, 4 ignored, 38 binaries.

PHASE 5J ITEM 1 IS DONE, AND NOT THE WAY IT WAS ATTEMPTED
  Four attempts went through the velocity field. The dominant transport at
  this scale needs no velocity field at all:

    D = kT / (6*pi*eta*r)        Stokes-Einstein

    20 nm bag, water, 320 K  ->  D = 2.1e-11 m2/s
    rms displacement in 1 s  ->  6.5 um = 6.5 cells of a micron grid

  Six cells per tick. Advection by a convection roll is a rounding error
  beside it. It is why a bacterium needs a flagellum: below a few microns,
  swimming loses to being shoved.

EVERYTHING IS DERIVED
  Temperature from the same derive() physics uses. Viscosity from the
  material table. Radius from the bag's own membrane area, r = sqrt(A/4pi),
  which it has because of the surfactant its chemistry made. A bag that
  grows slows down as 1/r, because Stokes said so.

THE LIMIT
  Beyond sigma ~ 1 cell/tick the walk saturates - a lattice cannot represent
  a bag crossing more than a cell in one hop. Stated at the code site.

WHAT IT UNBLOCKS
  Washout. A bag pinned to its birth cell could never reach a drain, so a
  population had births and no removal. Bags can now reach a Port::Sink,
  making the chemostat criterion testable: a lineage persists iff its
  division rate exceeds the dilution rate. That is where selection stops
  being enabled and starts being observed.

NOTE ON THE DETOUR
  The four failed attempts produced real fixes that stand on their own:
  implicit viscosity (conserving exactly at one sweep), and sor_omega being
  wrong by 100x on any grid that was not 41 cells. Both shipped separately.
