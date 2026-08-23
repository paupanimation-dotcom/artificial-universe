Artificial Universe — membrane transport stability limiter
==========================================================

4 files. Extract over your project and overwrite.

VERIFIED
  cargo test --release -p au-sim --test protocell -> 10 passed, 0 failed,
  1 ignored (see below). Full suite was 227 passed / 1 failed before the
  ignore was applied; that 1 is the ignored test.

THE DEFECT
  Membrane uptake was an explicit Euler step with no stability limit. A bag
  equilibrates in ~0.15 s against a 1 s tick, so every transfer overshot
  equilibrium ~7x. A daughter whose parent had eaten its substrate saw a
  near-infinite gradient and took 1563 molecules into an 837-molecule
  osmolyte budget, which burst it.

  before: 3000 nucleated, 2 divided, 2999 lysed, 3 alive
  after:  3000 nucleated, 1 divided, 0 lysed, 3001 alive

THE PRINCIPLE
  Transport may be wrong, but it may not be impossible. Fourth time this
  engine has needed it (Phase 2 CFL, 5b donor's pocket, 4b convex
  combination). A new transport path is not exempt for being new.

ALSO IN HERE: THE CLOCK
  physics.implicit_conduction was never set in any protocell world, so the
  physics layer capped the tick at 3.6e-4 s. Every protocell result so far
  was measured in a world about ONE SECOND old. Now enabled in both the
  scenario and the tests. 3000 ticks = 0.1 s wall clock, so 1e6-tick runs
  are cheap and we have been running far too short.

NOT VERIFIED - DO NOT RELY ON
  Protocell mobility (bags riding the fluid) is written and compiles, but
  its test is #[ignore]d: the scenario imposes a uniform current in a sealed
  tube, which incompressibility forbids. Needs periodic walls on the flow
  axis. bags_do_not_move_in_still_water also proves nothing right now - it
  passes when the velocity field is entirely absent.
