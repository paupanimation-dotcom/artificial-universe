# Changelog

## Protocells move — by Brownian motion, which needed no fluid solver at all

**250 tests passing, 0 failures, 4 ignored, 38 binaries** (was 246). Phase 5j's
first item, blocked for three sessions, is done — and the answer was not the one
being chased.

### Three sessions chasing the wrong physics

Protocell mobility was attempted four times through the velocity field. Each
attempt found something real and none produced a moving bag: a uniform current
cannot persist in a sealed tube; `physics.fluid.enabled` defaults false; a flow
cannot be declared, only driven; and a driven convection cell is too slow to run.
That last one sent this project into implicit viscosity, then into measuring the
pressure projection, then into finding `sor_omega` wrong by a hundredfold — all
real fixes, none of which moved a protocell.

**The dominant transport at this scale needs no field at all.** A vesicle is a
colloid, a colloid is kicked by the molecules around it, and Stokes–Einstein says
how hard:

```
    D = k_B·T / (6·π·η·r)

    20 nm bag, water, 320 K   ->   D = 2.1e-11 m²/s
    rms displacement in 1 s   ->   6.5 µm  =  6.5 cells of a micron grid
```

**Six cells per tick.** A convection roll is a rounding error beside it. This is
precisely why a bacterium needs a flagellum: below a few microns, swimming loses
to being shoved.

### Derived, not declared

Temperature comes from the same `derive` physics uses. Viscosity comes from the
material table. The radius comes from the bag's own membrane area — `r = √(A/4π)`
— which it has because of how much surfactant its chemistry made. **A bag that
grows slows down, as `1/r`, because that is what Stokes said.** Nothing in the
mobility path is a chosen number.

Same lattice treatment as advection: a hop probability per axis, so the mean is
right and no fractional cell offset accumulates as float state. Same derived RNG,
so a resumed world wanders identically — `brownian_wandering_is_reproducible`
pins the hash.

### The limit, stated

Beyond `σ ≈ 1` cell per tick the walk saturates: a bag crossing more than a cell
per step cannot be represented by a single hop. That is a resolution limit of a
lattice, it is where the model stops being quantitative, and it is written at the
code site rather than left to be discovered.

### What this unblocks

Washout. A bag pinned to its birth cell could never reach a drain, so a
population had births and no removal and could not turn over. Bags can now reach
a `Port::Sink`, which makes the chemostat criterion — a lineage persists iff its
division rate exceeds the dilution rate — testable for the first time. That is
the next thing, and it is where selection stops being enabled and starts being
observed.

### Files

* `au-sim/src/systems/protocell.rs` — Stokes–Einstein wandering, medium viscosity
* `au-sim/tests/protocell.rs` — +2 tests, a micron tube with no fluid solver

## The relaxation factor was a declared constant, and it was wrong by 100×

**au-physics: 51 passing, 0 failures, 1 ignored** (was 49). Every anchor still
holds, Ra_c included.

### The measurement that started it

Chasing why fine-grid flow is unaffordable, a probe found the pressure
projection is very nearly the *whole* cost of a fluid step — 20 ticks of a 24×12
convection cell take 4.09 s at sixty sweeps and 0.20 s at fifteen. Twentyfold,
for a fourfold cut.

The cheap win is not available, and measuring that mattered more than the
speedup would have. Divergence as a fraction of the velocity scale:

```
    15 sweeps   200.5%     grossly compressible
    30 sweeps    91.0%
    60 sweeps     1.4%     honest
```

An under-converged projection leaves matter accumulating where nothing put it,
and the buoyancy from that drifted density looks *exactly* like convection —
the artefact that gave Phase 2b fictitious stirring at half the critical
Rayleigh number with the perturbation set to zero. Sixty sweeps is near the edge
of honesty, not padding.

### And then the actual defect

Optimal over-relaxation for a Laplace problem is **not a constant**. It is
`2/(1+sin(π/n))` — a property of the *grid*, rising toward 2 as it refines:
1.59 at n=12, 1.77 at n=24, 1.86 at n=41, 1.91 at n=64.

`sor_omega` was declared 1.85, which is the optimum for a **41-cell** grid: the
tank Phase 2b validated Ra_c on. It has been silently wrong for every grid since.
Measured on a 12-cell grid at twenty sweeps:

```
    omega 1.00  ->  83.4% divergence
    omega 1.40  ->  15.4%
    omega 1.59  ->   1.8%     <- derived
    omega 1.85  -> 185.1%     <- the declared constant
    omega 1.95  -> 168.3%
```

**The declared value was worse than no over-relaxation at all.** Past the
optimum SOR does not converge slowly; it oscillates. A hundredfold error, in a
number nobody derived — which is precisely what §21b says not to do, written
before that section existed.

`sor_omega` now defaults to 0, meaning *derive it*, from the smallest active
dimension because convergence is set by the slowest direction to communicate
across. A nonzero value is still honoured for tests that need to pin one.

### What it did not buy, stated plainly

The convection cell got *slower*: 5.34 s against 4.09 s for twenty ticks. The
hypothesis is that an honest projection lets the flow actually convect, so
velocities rise and the advective Courant bound tightens — the earlier run being
fast because the flow was garbage. **That is a hypothesis and it is untested.**
Correctness improved and measured throughput did not.

### Where the speed actually has to come from

Single-grid Gauss–Seidel on a Laplace problem converges at a rate set by the
grid, roughly `cos²(π/n)` per sweep. That is the identical finding Phase 5b
recorded for deep-time diffusion, where the fix was `2n²` sweeps and a note that
multigrid is the real answer. It is the real answer here too, and one multigrid
would serve conduction, species transport and the projection together.

### Files

* `au-physics/src/grid.rs` — `derived_sor_omega`, `sor_omega` defaults to derive
* `au-physics/src/fluid.rs` — the derivation, and the cost measurements recorded
  at the projection itself
* `au-physics/tests/fluid.rs` — +2 measurement tests

## Implicit viscosity — built, conserving, and worth less than expected

**au-physics: 49 passing, 0 failures, 1 ignored** (was 47). Off by default
(`physics.viscous_iters = 0`), so every earlier phase runs byte-identically.

### What it is

Backward Euler on the viscous term, solved by Gauss–Seidel per velocity
component — the last member of the trio Phase 4b began, and the thing
`au-physics/src/implicit.rs` said at its head it was written to enable.

The integration is smaller than the roadmap assumed, because Phase 4b's
discipline carries over intact: **the solver proposes, the transport disposes.**
`viscous_solve` computes a better velocity field and moves nothing. The existing
symmetric-flux loop then runs *unchanged*, reading that field instead of the
current one. One word of difference.

Plain Gauss–Seidel, `omega = 1.0`, never over-relaxed: the update is a convex
combination and over-relaxation forfeits the maximum principle. Phase 4b paid
for that lesson with 4244 K in a box that started between 400 K and 1200 K.

### The property that makes it safe

`an_underconverged_viscous_solve_still_conserves_momentum_exactly` runs the solve
crippled to a **single sweep** for 500 steps and demands the books balance. They
do, exactly, because the two halves of a symmetric integer transfer are the same
integer.

Conservation does not depend on convergence. The trap this rules out is writing
the update as `Δp = m·(u_new − u_old)` per cell, which is conservative only on a
perfect solve and silently couples the engine's central invariant to an iteration
count.

### And the honest measurement

The point of building it was the protocell convection cell, where explicit
viscosity substepped for eleven minutes without finishing 460 ticks. With the
implicit solve on:

```
  20 ticks, 24x12 convection cell
    explicit    7.8 s
    implicit    4.6 s        1.7x
```

**1.7×, not orders of magnitude.** Viscosity was never the dominant cost there,
which the roadmap and I both assumed without measuring. At 0.23 s/tick the test
still needs about two minutes, so protocell mobility remains unverified for the
same performance reason as before.

What is left is structural: the advective Courant bound `dt ≤ dx/u`, which **no
implicit viscous scheme removes** because advection is hyperbolic rather than
diffusive, and the pressure projection's sixty sweeps per substep. Both should be
measured before anything is built on the assumption that fine-grid flow is now
affordable. It is cheaper; it is not cheap.

### Also unverified, and marked

`implicit_viscosity_takes_the_viscous_term_off_the_clock` is `#[ignore]`d because
**the test is wrong, not the code**: the Rayleigh–Bénard tank is calibrated to
Prandtl 1 at Ra≈2000, where a tick already fits inside the stable bound, so
explicit and implicit both take one substep and the comparison measures nothing.
The exemption is written and reaches the limiter; the claim that it removes the
limit has no case yet where viscosity binds.

### Files

* `au-physics/src/fluid.rs` — `viscous_solve`, the proposal step, scratch buffers
* `au-physics/src/grid.rs` — `Boundary::viscous_iters` (default 0)
* `au-physics/src/transport.rs` — viscous bound exempted when implicit
* `au-physics/tests/fluid.rs` — +3 tests
* `au-sim/src/systems/physics.rs` — `physics.viscous_iters`
* `au-sim/tests/protocell.rs` — enabled in the convection scenario; measurements

## Nucleation bounded by material, and a recurring error finally named

**246 tests passing, 1 ignored, 0 failures** (was 244). Two changes and one
piece of bookkeeping that is worth more than either.

### Nucleation was bounded by a `for` loop, not by physics

The protocell system closed *one* bag per cell per tick. That is not a rate
anyone derived — it is the shape of the loop — and a trace found the consequence:
the medium held 42,000,000 assembled surfactant against a bag's 12,567, enough
for 3,300 vesicles, drawn down one per tick until the pool emptied around tick
4,200. Bags therefore appeared long after the substrate they needed was gone.

The first diagnosis was wrong and reading the code disproved it. I assumed the
42M was a supersaturated solution that Phase 5f had failed to limit;
`assembled_count` caps free monomer at the CMC and treats the rest as surface,
so those molecules were *already structure* — micelle and sheet, not solution.
The model was right.

A sheet closes because a closed vesicle has no exposed edge, and nothing in that
says one per second. So every bag the surface can support now closes:

```
                    nucleated   divided   lysed
  one per tick           4130         1       0
  material-bounded       4130        64       0
```

`MAX_NUCLEATIONS_PER_CELL_TICK = 64` is **declared**, and it bounds *work* rather
than physics: surface left over closes on the next tick, so the only thing it
changes is how finely the work is spread.

### The new ceiling is competition, which is the honest result

64 divisions, identical across a tenfold range of starting surfactant — because
it is the *first tick's* bags that divide, before the medium fills with
thousands of competitors. After ~150 ticks the population is frozen for the
remaining 2,850: 4,000 bags sharing one cell's substrate, none with the
throughput to grow.

That is what a chemostat exists to break, and it makes **washout the binding
constraint** rather than a nice-to-have. Which needs mobility, which is below.

### Protocell mobility: three attempts, still unverified

Written, compiled, never demonstrated. Each attempt found a real thing and none
produced a moving bag:

1. A uniform current in a sealed tube — incompressibility forbids it, and the
   pressure projection correctly annihilated it.
2. `physics.fluid.enabled` defaults to **false**, so the physics system was
   zeroing momentum columns it did not own. The column existing is not the
   solver running.
3. The physics system gathers, solves and scatters, so its own state overwrites
   an initial condition. **A flow cannot be declared; it has to be driven** — by
   buoyancy or a pressure boundary, the way `au convect` has since Phase 2b.

Attempt four rebuilt the test as a real convection cell and hit the wall Phase 2b
predicted: a 24×12 fluid solve substepping to the viscous limit exceeded eleven
minutes without finishing 460 ticks. Both convection tests are `#[ignore]`d with
that written at the test. **Implicit viscosity is the blocker**, it is the last
member of the trio Phase 4b began, and `au-physics/src/implicit.rs` says at its
own head that it was written to be reusable for exactly this.

`bags_do_not_move_in_still_water` had been passing while the velocity field was
entirely absent. It now carries the control it needed, commented out alongside
the ignored tests so that it stops being decoration in one edit.

### A number means nothing without its container

The architecture doc gains a section naming a class of error this project has
hit **seven times since Phase 3** and documented as seven separate traps:

```
  2e9 molecules      → a µJ-quantised cell    → combustion released no heat
  a driving ΔT       → a coarser mass quantum → fictitious convection
  1e6 molecules      → the heating tests      → a test measuring a units bug
  a C12 alcohol      → a phase validated on C5→ seconds per tick
  a 1 s tick         → a 1 µm cell            → diffusion number 1000
  60,000 surfactant  → a cell 1e6× larger     → below the CMC, nothing assembles
```

Documenting the instances and never the pattern is why each new one looked
novel. The check is arithmetic and belongs before the literal, not after the
failure: a population against `CMC × volume` and the ~2×10¹² event floor; a
timestep against `D·dt/dx²`; a molecule against the size the phase below
validated at.

### Files

* `au-sim/src/systems/protocell.rs` — material-bounded nucleation
* `au-sim/tests/protocell.rs` — convection-driven drift tests, honestly marked
* `ARCHITECTURE.md` — "A number means nothing without its container"

## Membrane transport had no stability limit — and it was killing every protocell

**Every bag in the world was dying of a numerical error in the uptake step, not
of physics.** With the fix, mass lysis disappears entirely:

```
                    nucleated   divided   lysed   alive
  before                 3000         2    2999       3
  after                  3000         1       0    3001
```

### The defect

Uptake across a membrane was an explicit Euler step with no limiter. A bag
equilibrates with its medium in about 0.15 s; the tick is 1 s. Every transfer
therefore overshot equilibrium roughly sevenfold and reversed on the next tick.

The trace that found it:

```
UP  sp=pentane  c_in=0.0  c_out=2.0e25  flux=1563  n=1563
    sp=water    c_in=1.47e25  c_out=1.00e25  n=-173
    sp=H2       c_in=5.53e25  c_out=4.00e25  n=-649
```

A daughter inherits almost no substrate, because its parent's own chemistry
consumed it. That leaves a near-infinite inward gradient, and 1,563 molecules
arrive against an osmolyte budget of 837. The osmotic volume explodes and the
bag bursts. Reported last session as a finding about protocell fragility; it was
mostly this.

### The fix, and the principle it belongs to

Transfers are now capped at the amount that brings both sides to equal
concentration — `Δc · V_in · V_out / (V_in + V_out)`, which for a bag against its
cell is `Δc · V_in` to a part in a million.

This engine has now learned the same lesson four times: Phase 2's CFL
substepping, Phase 5b's donor's-pocket limiter, Phase 4b's convex combination,
and this. **Transport may be wrong, but it may not be impossible.** Every
transport path needs a limiter, and a new one is not exempt for being new.

### Why nothing caught it earlier

The clock. `physics.implicit_conduction` was never set in any protocell world, so
`PHYSICS_TIME_CAPPED` correctly refused the requested one-second tick and took
3.6e-4 s instead — which made the uptake step accidentally stable, and meant
**every protocell result to date was measured in a world about one second old**.
Both the scenario and the tests now enable the implicit solver. 3000 ticks runs
in 0.1 s of wall clock, so runs of 1e6 ticks are affordable and everything so far
has been looking through the wrong end of the telescope.

### Honest status

* **Bags are stable and inert.** Lysis is gone; growth is not there. Nucleation
  still floods at one per tick. That is now the real problem rather than an
  artifact sitting on two others.
* **Protocell mobility is written, compiled and UNVERIFIED.** Its test is marked
  `#[ignore]` with the reason: the scenario imposes a uniform current in a sealed
  tube, which incompressibility forbids, so the projection annihilates it and
  nothing drifts. Needs periodic walls. Do not rely on mobility.
* `bags_do_not_move_in_still_water` passes identically when the velocity field is
  absent, so it currently proves nothing. Labelled as such at the test. That is
  the third test in this project to pass because of a defect.

### Files

* `au-sim/src/systems/protocell.rs` — the equalisation limiter; bag advection
* `au-sim/tests/protocell.rs` — implicit conduction; drift tests; honesty labels
* `au-headless/src/main.rs` — implicit conduction in the `au watch` scenario

## Phase 5j groundwork — a drain for individuals, and what the view found

**244 tests passing, 0 failures, 38 binaries** (was 237). This is a small change
and a large finding; the finding is the point.

### What `au watch` showed on its first sweep

Six runs, varying the medium's starting surfactant and the substrate feed:

```
feed on,  surfactant 2.0e6 / 8.0e5 / 6.5e5  →  3000 nucleated, 2 divided, 531 lysed
feed off, surfactant 2.0e6 / 8.0e5 / 6.5e5  →  ~1950 nucleated, 0 divided, ~330 lysed
```

Two things, neither of which any test was going to report.

**The starting composition is irrelevant.** Every fed run gives byte-identical
counts across a threefold range of initial surfactant, because the bulk
manufactures its own from the fed substrate. The medium is a vesicle factory, and
what it starts with does not matter.

**Turn the feed off and nothing divides at all.** Bags need throughput; without it
they nucleate, sit, and burst.

So Phase 5i step 2's headline result — composition determines replication rate —
rests on two division events against zero. The claim is true and the magnitude is
tiny, which the changelog for that phase did not say because nothing could see it.

### The structural gap this exposed

Bags could only ever *die by bursting*. There was no way for one to leave. A
population with births and no removal cannot turn over, and without turnover a
faster-dividing lineage never displaces anything — it merely adds to a pile.
Selection needs death that is not failure.

Phase 5h gave matter a drain. Individuals never got one. `Port::Sink` now removes
protocells in its cell along with everything else, booking their contents as
outflow so the atom ledger still balances exactly. The drain stays blind: it
removes bags without regard to composition, which is what keeps selection an
outcome rather than an input.

### And why that is not enough, stated plainly

**The mechanism cannot currently fire on its own.** A `Port::Sink` holds its cell
at zero, so surfactant never accumulates there, so no bag ever nucleates in the
one cell that would remove it. Washout needs protocells to *move*, and they have
no mobility at all — a bag is pinned to the cell it formed in, forever.

`a_bag_in_an_open_port_washes_out_and_the_books_record_it` places a bag by hand
and proves the bookkeeping is right for the day mobility arrives. It is not
evidence that populations turn over, and it is labelled so at the test.

**Protocell mobility is therefore the first item of Phase 5j**, ahead of anything
about competition or selection. Bags should ride the fluid the way species ride
their gradients — Phase 2b's velocity field already exists and species transport
has never been coupled to it either, so the two are one piece of work.

### Files

* `au-sim/src/systems/protocell.rs` — sink cells, washout, gross booking
* `au-sim/tests/protocell.rs` — +1 test, with its own limits recorded
* `au-headless/src/main.rs` — `--surfactant`, `--pentane`, `--feed` for `au watch`

## Observation — a way to look at the world (and the first thing it found)

**`au watch` writes a self-contained HTML page of what a run did.** No
dependencies, no JavaScript, no CDN: inline SVG in one file that opens from a
`file://` URL. 243 tests passing (was 237); +6 for the recorder.

### Why this came before the next phase

Every result in this project has come from a test, and a test can only check
what somebody already thought to ask. That has worked for correctness — the
suite has caught a ledger netting opposing flows, a vesicle dissolving itself, a
reaction creating a carbon atom — and it is close to useless for discovery.

Three separate debugging sessions ended in an `eprintln!` in a hot loop and a
`grep`. Each was solved the moment real numbers appeared: a world cycling every
two ticks, bags whose contents never changed, a cell at 10¹⁸ K and then at
absolute zero, compartments starved by their own medium. Four time series would
have shown all four in seconds. The A-life literature is unanimous that finding
interesting behaviour requires watching the thing, and this project had no way
to watch anything.

### Read-only, enforced by the compiler

`Recorder::observe` takes `&World`, not `&mut World`. It cannot write a column,
emit an event, touch the allocator, or change the hash — not by convention but
because the borrow checker refuses. It lives outside the scheduler; no system
runs it, and nothing in the layer stack depends on it.

`observing_a_world_never_changes_it` runs two identical worlds, observes one
twice per tick for three hundred ticks, and asserts both the hash and the saved
bytes match. Phase 1 shipped a bug where an event emitted during chunk eviction
made *RAM pressure* part of world identity; that is the failure this rules out.

### What it draws, and why those four

Each panel corresponds to a mistake that actually happened unseen:

* **the protocell census** — living population against cumulative births,
  divisions and deaths. A sawtooth with rising counters is churn, not growth,
  and neither series shows that alone.
* **populations, medium against encapsulated** — drawn apart, never summed.
  Since protocells, matter lives in two places; their sum is what conservation
  checks and their *ratio* is the only thing that says whether compartments are
  doing anything. Log scale, because these span orders of magnitude.
* **temperature as a band** — on a signed log scale, so 10¹⁸ K stays on the page
  instead of flattening everything else into the axis.
* **matter across the boundary** — gross, never net. Equal and rising is a
  steady state; equal and flat is a sealed box; diverging is accumulation.

### Instrumentation is not rendering

Stated in the module docs because it would be easy to blur. Rendering proper
must **derive** appearance: a colour will be a consequence of what a thing is
made of, and the human view will be one declared mapping among possible sensors.
Nothing here does that. These are arbitrary colours assigned to arbitrary series
so a person can tell two lines apart, drawing *measurements of* the world rather
than the world. The page says so at the bottom.

### The first thing it found, in one run

```
3000 ticks · 2471 protocells · 3000 nucleated · 531 lysed · 2 divided
```

**One nucleation every tick, forever, and two divisions in the whole run.**

The population is not a population of reproducing bags. It is a conveyor belt of
freshly nucleated ones, most of which never live long enough to grow, plus a
steady trickle of lysis. Phase 5i step 2's headline test compares productive
against inert division counts and passes on `2 > 0` — the claim it makes is true,
and the *magnitude* behind it is far weaker than the changelog entry implied.

That is the honest correction, and no test in the suite was going to produce it.
Nucleation is unbounded because the medium is continuously resupplied — by the
pentane port and by bulk chemistry making more surfactant — so a bag that would
divide is drowned by newborns. Recorded as the first item for Phase 5j, which is
where populations were always going to become the subject.

### Files

* `au-sim/src/observe.rs` — `Recorder`, `Sample`, `BagSample`; read-only
* `au-sim/tests/observe.rs` — 6 tests, including the no-perturbation guarantee
* `au-headless/src/render.rs` — HTML/SVG, zero dependencies
* `au-headless/src/main.rs` — `au watch`, and the built-in protocell scenario

## Units fix — reaction enthalpy was molar, spent per event

**237 tests passing, 0 failures, 37 binaries.** Same count as before: a defect
was fixed and the two tests that were measuring it were corrected.

### The defect

`BondEnergyModel::base_per_order` is quoted in **J·mol⁻¹**, because that is how
bond energies are tabulated. Everything built from it is molar. But
`reaction_enthalpy` promised *per event* in its doc comment, and two of its three
callers believed the prose:

```
  network.rs        dh_molar → / AVOGADRO            correct all along
  main.rs (demo)    .as_joules() / AVOGADRO          correct all along
  reactor.rs        rxn.enthalpy_j = e.as_joules()   molar spent per event
  au-sim/chemistry  enthalpy_j = ...as_joules()      molar spent per event
```

Every declared reaction without an explicit `enthalpy_j` delivered **6.022×10²³
times too much heat**. It is what drove a 10⁻¹⁵ kg cell to absolute zero in
Phase 5i step 2, where it was worked around by declaring `enthalpy_j` and logged
rather than chased.

### The fix

The unit now lives in the name and in the type. `reaction_enthalpy_molar`
returns a bare `f64` of joules per mole rather than a quantised `Energy`, so the
conversion cannot be skipped by accident, and both wrong callers now divide.
`network.rs` is unchanged — it was right from the start, and named its variable
`dh_molar` to say so.

### The lazy fix, and the anchor that refused it

Dividing inside `reaction_enthalpy` itself was tried first. It broke
`generated_reactions_reach_boltzmann_equilibrium` immediately — measured K
collapsed to exactly zero against an expected 8.9554e-29 — because `network.rs`
would then have divided twice. That expected value implies a per-event ΔH of
~1.78e-18 J, so the anchor already assumed the conversion happened downstream.
The test was right and the patch was wrong; it was reverted, documented at the
line it concerned, and the call graph traced properly instead.

### Two tests were passing *because* of the bug

`Energy` is quantised to the microjoule and one bond is ~5×10⁻¹⁹ J, so a
**correct** per-event enthalpy rounds to exactly zero in an `Energy`. Two
assertions read `Reaction::enthalpy` and only worked because a molar value is
large enough to survive quantisation. One of them legitimately builds its
reactions by hand at molar scale and now says so; the other goes through the
derived path and reads `enthalpy_j`.

The heating tests needed real populations. Below roughly **2×10¹² reaction
events** a microjoule-quantised thermal field cannot resolve molecule-scale
chemistry at all — one microjoule divided by one bond. The au-chem test went from
10⁶ molecules to 10¹⁸ and the sim-loop test from 4×10⁶ to 4×10¹³. That is
arithmetic, not a fudge: the coupling is real and simply needs a real number of
molecules.

### Why six phases missed it

**Energy conservation held perfectly the entire time.** The heat was exactly
accounted for on both sides of the ledger — it was the wrong size. A conservation
check cannot detect a dimensional error, because the error is consistent
everywhere it appears.

This is the same shape as the Phase 5h ledger that netted opposing flows and
reported a sealed world in the middle of a torrent: the invariant was true and
the world was still wrong. **Invariants catch bookkeeping mistakes; only an
external anchor catches a dimensional one.** It also survived because most
declared-reaction tests supply `enthalpy_j` explicitly and never reach the line.

### Files

* `au-chem/src/reaction.rs` — `reaction_enthalpy_molar`, returns `f64` J·mol⁻¹;
  `Reaction::enthalpy` documented as unusable for bond-derived reactions
* `au-chem/src/reactor.rs` — `with_derived_enthalpy` converts to per event
* `au-chem/src/network.rs` — call-site rename only; behaviour unchanged
* `au-sim/src/chemistry.rs` — declared-reaction path converts to per event
* `au-headless/src/main.rs` — demo reads `enthalpy_j`, scales up for kJ/mol
* `au-chem/tests/reactions.rs`, `au-sim/tests/chemistry.rs` — corrected, with
  the reason recorded in each

## Investigation — a units defect in the declared-reaction energy path (no behaviour change)

**One file, comments only. 237 tests still pass; nothing is altered at runtime.**
This entry exists because the finding is worth more than the fix would have been,
and a half-understood correction to the energy model is worse than a documented
defect.

### What was found

`BondEnergyModel::base_per_order` is documented in **J·mol⁻¹** — a C–C single
bond is ~350 kJ/mol — so `molecule_energy` returns a molar quantity.
`reaction_enthalpy` subtracts two of those and hands the result straight to
`Energy::from_joules`, while its own doc comment promises a value **per event**,
because what it returns is added to one cell's thermal energy each time one
reaction fires.

Those two cannot both be true. `network.rs` has divided by Avogadro's number
since Phase 5a, so the **generated** path is correct; the **declared** path
appears to be wrong by 6.022e23.

This surfaced from Phase 5i step 2, where a declared reaction
`C₅H₁₂ + H₂O → C₅H₁₂O + H₂` drove a 10⁻¹⁵ kg cell to absolute zero. The test
worked around it by declaring `enthalpy_j`; this is the follow-up that should
have happened before anything else was built on top.

### Why it was not fixed

Dividing by Avogadro here breaks `generated_reactions_reach_boltzmann_equilibrium`
— measured K collapses to exactly zero against an expected 8.9554e-29. That
expected value implies a per-event ΔH of ~1.78e-18 J, so **the anchor test
already assumes the conversion has happened.** Some consumer is therefore reading
this function's molar output and converting downstream, and dividing at the
source converts twice.

The correct fix is to find that consumer and make the units explicit end to end.
That is careful work on the layer everything else stands on, and it is not
something to rush. Reverted deliberately, with the reasoning written at the line
it concerns so the next attempt starts from the evidence rather than repeating it.

### Why six phases missed it

**Energy conservation held the entire time.** The heat was exactly accounted for
on both sides of the ledger — it was simply the wrong size. A conservation check
cannot detect a dimensional error, because the error is consistent everywhere it
appears.

This is the same shape of failure as the Phase 5h ledger that netted opposing
flows and reported a sealed world in the middle of a torrent: the invariant was
true and the world was still wrong. **Invariants catch bookkeeping mistakes;
only an external anchor catches dimensional ones.** The anchor here is arithmetic
anyone can do by hand — rearranging a few covalent bonds moves on the order of
1e-19 J, not 1e5.

It also survived because most declared-reaction tests supply `enthalpy_j`
explicitly and never reach this line.

### Standing item

Logged against Phase 3b. Until it is resolved, declared reactions should state
`enthalpy_j` rather than relying on the bond-derived value.

## Phase 5i (step 2) — Protocells: composition determines replication rate

Bags are now **individuals**. They nucleate out of the medium, feed through
their own membranes, run their own chemistry at their own concentration, divide
when geometry allows, and record which bag they came from.

### The claim, and why it is the only one that mattered

A compartment that merely persists is a rock. Phase 5h showed a membrane is
worth having in a flow — 457× retention against a blind drain — and step 1 showed
when a bag has enough surface to become two. Neither gives anything to inherit.
Without a link from *what a bag is made of* to *how fast it divides*, composition
is inherited information with no consequences, variation cannot affect
persistence, and no adaptation is possible. That is the crystal-with-a-membrane
failure, and it is the one artificial life hits most often.

The link is now closed, and the trace is the result:

```
events=335   inner0=335   surfactant inside: 12,567 → 12,902
BAG v=0.6605  area=5.295e-15      ← below 1/√2, fission
BAG v=0.9030  area=2.676e-15      ← daughter: area halved, a stable sphere again
```

A bag imports pentane, converts it to surfactant *inside*, and the surfactant
cannot leave — because it **is** the membrane. Area ratchets upward, `v` crosses
the two-sphere bound, the bag divides, and both daughters inherit the composition
that caused it. A world whose chemistry cannot build membrane never divides at
all. Nothing anywhere grants an advantage; the advantage is what trapping a
product does in a compartment.

The daughter landing at v ≈ 0.90 is the √2 jump step 1 predicted from geometry,
observed rather than asserted.

### Five bugs, and what each one cost

**Nucleation built one bag out of the whole medium.** With 40M assembled
surfactant, a single vesicle of radius 1.38 µm formed inside a 1 µm cell — eleven
times its container's volume. It swallowed everything, fissioned, and its
daughters burst because the emptied medium made their osmotic volume explode; the
contents dumped back and it renucleated, forever. A two-tick cycle that no
chemistry could influence, which is exactly why both test worlds divided 300 times
each and looked identical. The error was physical: surfactant above the CMC forms
**many small vesicles**, not one large one, because bending energy is what would
pay for a large one. Nucleation now takes only the molecules needed to close the
smallest bendable bag and leaves the rest. Bags are born as spheres at `v ≈ 1` and
can only reach fission by growing.

**A vesicle dissolved itself.** Membrane species were allowed to diffuse out
through their own membrane into a medium sitting at the CMC. Phase 5g had already
made this exclusion for grid transport — structure has no independent gradient to
run down. The consequence is deliberate: a bag cannot grow by absorbing monomer,
so the only route is making surfactant itself.

**A nearly-empty bag hung the reactor.** Near-zero osmotic volume means contents
look infinitely concentrated, and the reactor correctly subdivides the timestep
toward infinity to hold its consumption cap. Fixed with an excluded-volume floor:
molecules occupy space, so a bag of N of them is at least `N·v_molecule` across.

**The save frame could not be read backwards.** The life trailer's lengths were
written at the *front* of its section, and the reader then indexed from the front
of the remaining bytes — which is the registry, not the store. Every length now
sits at the end, so the section can be located walking back from the magic. The
atom trailer is written first and the life trailer last, so `load` peels life,
then atoms, then finds the registry.

**Two test chemistries the engine was right to reject.** Atomic oxygen is a
diradical with two unsatisfied valences; forming bonds from it is enormously
exothermic, and five million such reactions heated a 10⁻¹⁵ kg cell to 10¹⁸ K.
Water balances the same reaction and is closed-shell — and then the bond-derived
enthalpy drove the cell to absolute zero instead. **Twenty million molecules
reacting in a cubic micron is a macroscopic energy event, and a model with exact
thermodynamics says so rather than quietly averaging it away.** The reaction
enthalpy is now declared rather than derived, for the reason below.

### Open question logged against Phase 3b

Left to the bond model, `C₅H₁₂ + H₂O → C₅H₁₂O + H₂` produces a per-event enthalpy
roughly **eleven orders of magnitude** larger than the bond balance predicts
(hand-computing C–C, C–H, O–H, C–O and H–H gives about +75 kJ/mol, or 1.2e-19 J
per event; the observed total implies ~8.7e-8 J per event). Whether the fault is
in the bond model or in these particular graphs is not a protocell question, and
a test of compartmentalisation should not be silently measuring it. `enthalpy_j`
is declared in the test and the discrepancy is logged.

### The medium outcompetes the compartments, and the open system fixes it

Bulk and bags run the same chemistry and the bulk has a million times the
material, so compartmentalisation buys no advantage in *rate* — only in trapping
the product. The medium converted every substrate molecule before a single
vesicle formed, and the bags nucleated into an exhausted world. A Dirichlet
source from Phase 5h holds the substrate at fixed concentration so consumption
cannot exhaust it. Both worlds run under identical conditions; the only
difference remains whether chemistry inside a bag can build membrane.

### Structure

* `Protocell` — id, parent, birth tick, chunk, cell, sparse contents. Division
  frees the parent's id and allocates two fresh ones, each recording the parent,
  so the lineage is a clean binary tree and no id ever means two things. This is
  the first use of the `EntityId` generation counter that Phase 1 added with the
  note that a family tree must never point at a stranger.
* `ProtocellStore` — flat, sorted by id, so iteration order is the same
  everywhere and forever.
* **Species populations now live in two places.** Grid columns hold the bulk;
  protocells hold what they have taken up. Every accounting must sum both, and
  `matter_is_conserved_across_the_membrane` checks it every tick.

### Honest limits

* **Declared reactions only inside bags.** Open-ended generation needs the
  reaction cache shared between two systems rather than rebuilt in each.
* **Reaction heat is dumped into the chunk's first cell.** Exact, but crude.
* **Selection is enabled, not observed.** One bag dividing faster than another is
  the precondition. Watching a composition *spread through a population* is
  Phase 5j.
* **`clock.rescale.enabled = false` did not stop rescaling** — `dt` came out at
  3.6e-4 rather than 1.0. Harmless here, unexplained, logged.

## Phase 5i (step 1) — Division, and heredity without a genome (incremental update)

**230 tests passing** (was 215). Nothing is wired into the simulation yet: this
step is pure `au-chem`, and every world runs byte-identically to Phase 5h. What
it establishes is *when a bag divides* and *what its daughters get* — the two
questions that have to be answered before protocells can exist as objects.

Phase 5h left compartments with a reason to exist: in a flow, a barrier keeps
its contents, by a factor of 457. But a compartment that merely persists is a
rock. It has no descendants, so nothing about it can be inherited, and the
difference between a good barrier and a bad one cannot accumulate.

### Nobody chooses a size threshold

The tempting implementation is a rule — divide when the contents exceed N — and
that number would instantly be the most load-bearing constant in the engine.
The Golden Rule forbids exactly this. So the criterion is taken from geometry,
where it has been sitting all along.

A bag has an **area** set by how many amphiphiles it has assembled (`A = N·a`,
Phase 5f's constant, unchanged) and a **volume** set by osmosis: a membrane
passes water and not solute, so the bag swells until the concentration inside
matches the concentration outside, `V = n_internal / c_external`. Not a knob —
a consequence of what the bag holds and what surrounds it.

Their ratio is the reduced volume, `v = 3V√(4π)/A^(3/2)`, the standard order
parameter of vesicle shape (Seifert, Berndl & Lipowsky, 1991). Two facts about
it decide everything, and both are arithmetic:

* **v > 1 is impossible.** A sphere encloses the most volume per unit area
  there is, so a bag whose osmotic volume exceeds its spherical maximum cannot
  hold itself together. Real bilayers survive about 3% areal strain first and
  then fail.
* **v ≤ 1/√2 admits two spheres.** Two spheres of radius r have area `8πr²` and
  volume `(8/3)πr³`; substitute and the reduced volume is exactly `1/√2 ≈
  0.7071`. Above it, symmetric fission needs more membrane than the bag owns.

`1/√2` is not tuned. It is what `√π/√(2π)` equals.

**The consequence nobody wrote down:** a protocell now sits between two
failures. Make osmolytes faster than membrane and it bursts; make membrane
faster than osmolytes and it divides. Whether a bag divides, bursts or persists
is decided by **the ratio of two rates in its own chemistry** — not by its size,
not by a timer. `fate_is_decided_by_a_ratio_of_two_chemistries_not_by_size`
checks the same three fates at 200,000, 1,800,000 and 40,000,000 amphiphiles;
`two_bags_of_very_different_sizes_share_a_fate_if_they_share_a_ratio` closes it.

### The determinism decision, made deliberately

Partition is the engine's **first stochastic physics**, and it is the first time
the question has been unavoidable. Every process so far had a mean-field answer.
Partition does not: when a bag pinches in two, which side a molecule ends up on
is decided by where it happened to be, and there is no deterministic function of
the parent that gives it. Phase 5b met the same question at the diffusion
quantisation floor and deferred it explicitly. This is where it gets answered.

* **Variation has to come from somewhere.** Evolution is variation plus
  differential persistence. 5h supplied the second. Without the first, every
  daughter is its parent forever and nothing can happen.
* **Determinism survives untouched.** The engine has no global RNG: randomness
  is derived, `f(seed, domain, where, when)`, so a resumed world draws identical
  numbers. What arrives is *stochasticity*, not nondeterminism. The two are
  routinely confused and are not the same thing.

The draw is exact below 1024 copies — for a fair coin, a population count of
random bits, which is what a binomial is — and normal-approximated above, where
the error is already under `1/√n` and the quantity being approximated is itself
the noise. The switch is stated, not hidden.

### Heredity without a genome, measured

No template. No copying. No sequence. Only a bag splitting in two — and a
daughter resembles its parent **fifty times more closely** than it resembles a
composition built from the same abundances rearranged. The proportions survive
because the coin is fair *per molecule*. This is compositional inheritance
(Segré, Ben-Eli & Lancet, 2000): heredity carried by what a thing is made of
rather than by anything written down.

And fidelity is not a setting. A species present in `n` copies is transmitted
with relative error `1/(2√n)`, so:

```
  same recipe, four scales        mean compositional distance parent → daughter

     ×10                                    (worst)
     ×1,000                                    ↓
     ×100,000                                  ↓
     ×10,000,000                          < 1e-8
```

**Copy number is heredity's fidelity**, and it falls out of counting rather than
being declared — which is the same reason a real cell carries one genome and not
one molecule of each protein.

### The honest limit, written as a test

`compositional_inheritance_decays_across_generations` exists to stop this being
oversold. Take a daughter, grow it back, divide again, repeat: distance from the
founder grows every generation, because each split adds its own √n noise and
nothing corrects it. There is no proofreading — a composition cannot be
*checked*, only re-drawn. This is the known weakness of compositional genomes
(Vasas, Szathmáry & Santos, 2010): they carry information but cannot hold it
still, so selection has very little to grip.

Recording it now means step 2 knows what problem it is solving instead of
discovering it in three phases' time. A template — something that can be
compared against and corrected — is what fixes this, and that is Phase 6.

### Also honest

* `1/√2` is the *geometric* bound: where two daughters become possible, not
  where they become energetically favourable. The real budding transition
  depends on bending rigidity and spontaneous curvature and sits higher. The
  energetics would move the line, not remove it.
* **Nothing is wired in.** No protocell exists as an object yet; there is no
  storage, no identity, no lineage. Step 2 is entities, and it is the larger
  half: variable-size composition per protocell needs a codec extension, and
  `EntityAllocator` — built in Phase 1 with the comment that a family tree must
  never point at a stranger — finally gets used for what it was written for.

### New surface

* `au-chem/src/vesicle.rs` — `shape`, `reduced_volume`, `fate`, `Fate`, `Shape`,
  `partition`, `composition_distance`, `FISSION_REDUCED_VOLUME`
* Tests: +15 (`au-chem/tests/vesicle.rs`)

## Phase 5h — The open system: matter can cross a boundary (incremental update)

**215 tests passing** (was 196). Every earlier test passes unchanged; a world
without ports (the default) runs byte-identically to Phase 5g.

Energy has been able to cross a boundary since Phase 2 — heat in at the floor,
radiation out at the top, and everything interesting living in the gap. Matter
never could. The reactor, the vent and the RAF survey all took a fixed inventory
of atoms and rearranged it, which is why Phase 5c was careful to claim only that
its autocatalytic set *could* sustain itself **if it were fed**. Nothing had ever
fed anything.

A box has one destiny: it runs down to equilibrium and stops. Building
individuality and division on top of that would have meant tuning the rules of
growth against a chemistry that was dying anyway, then retuning them once flux
arrived. So the ordering was: open the system first.

### Dirichlet, not injection

The obvious implementation is a source term — add N atoms per second — and it is
wrong twice. It needs a fractional accumulator to avoid a deadband at low rates
(a species whose expected feed is 0.3 molecules per tick would otherwise never
arrive at all), and it **declares the flux**, when the flux is precisely the
thing that ought to be derived.

So a port is a Dirichlet condition on population, the direct twin of the
`bottom_temp` / `top_temp` conditions the fluid solver has used since 2b:

* `Port::Source` — a window onto a reservoir held at a fixed composition, which
  tops its cell **up** and never down;
* `Port::Sink` — a window onto an infinitely dilute exterior, which draws its
  cell to zero.

Nothing declares a rate. The port states a concentration, the gradient does the
rest, and the throughput is an outcome of geometry and diffusivity. Every
operation is an exact integer clamp: no rounding, no residue, no new state to
snapshot, and applying a port twice in a row is a no-op the second time.

**A sink may not have taste.** It removes every species without preference. A
drain that spared some molecules and not others would be a fitness function
written into a boundary condition — selection supplied as an input rather than
obtained as a result. This is the single most important constraint in the phase.

### The result that made the phase worth doing

Since Phase 5g, assembled surfactant lowers a cell's permeability, and
permeability multiplies diffusive exchange. Since this phase, a drain removes
whatever reaches it. Neither of those facts is about advantage. Put them
together and one appears anyway:

```
au-sim/tests/membrane.rs — same tube, same tracer, same 20,000 ticks,
                           an open port at one end

  tracer surviving the drain     with a membrane   328,945
                                 without             720
                                                   ─────────
                                                     457×
```

Nothing in the engine grants membranes an advantage. The advantage is what a
barrier *is*, in the presence of a flow — implied by two phases separately and
invisible until matter was allowed to leave.

This is **not yet selection**: there is no heredity, no individuals, nothing that
can differentially reproduce. It is the precondition — an environment in which
two chemistries of identical composition have measurably different prospects,
decided by physics rather than by a parameter. Phase 5i is where something can
inherit that difference.

### The books, and the bug that proved they needed rewriting

`AtomLedger` counts what has crossed, **per element** rather than in bulk mass: a
bug that turned carbon into oxygen would balance to the attogram and still be a
catastrophe. The identity, checked every tick in `the_books_balance_exactly_per_element`:

```text
    atoms_of(e, now) - atoms_of(e, start)  ==  in[e] - out[e]      EXACTLY
```

Conservation does not weaken when the world opens. "Totals never change" was
never the real invariant — it was the special case of this one where nothing
crossed.

**The bug.** The first ledger accumulated one *signed* delta per species and
booked that. At steady state a source feeding 2,490 molecules and a sink draining
2,490 sum to zero, so the books recorded a sealed world while matter poured
through it at full rate. What makes this worth writing down is that the net
identity above **still held** — both sides were zero — so the invariant that
everything else is checked against could not detect it. The only symptom was a
reported throughput of nothing in a tube that was visibly flowing. Gross in,
gross out, always; `opposing_ports_are_booked_gross_not_net` pins it.

**The second bug**, caught by an existing test rather than a new one: a sealed
world writes no atom trailer, so a resumed world got a zero-length ledger where a
fresh one had a zeroed ledger sized to the element table — different lengths,
different hash, and `catalysed_worlds_are_reproducible_and_resumable` failed on
the spot. An empty ledger and a zeroed ledger are the same statement about
history and must be the same world.

### Validation

The anchor is Fick's steady solution in one dimension, the twin of the law Phase 2
validated conduction against: **steady diffusion between a fixed source and a
fixed sink gives a linear profile**, with the same flux across every face. A
21-cell tube run to steady state lands on the analytic ramp to within 2% of the
held value, holds its contents to within 1% while continuing to pass matter
through, and returns what it is fed to within 2% face by face.

The distinction those three assertions draw together is the whole point:
**equilibrium is where a box ends; a steady state is somewhere a system can
live.**

### Persistence

The ledger is history, so it must survive a save; the ports are a rule, so they
live in config exactly as chemistry and materials do. The ledger is appended as a
trailer identified by a magic at the very end of the save frame — read backwards,
costing one comparison, leaving the existing layout untouched, and letting a save
written before this phase load with an empty ledger, which is the correct reading
of a world where nothing could cross.

### Honest limits

* **Ports are addressed by cell index within every active chunk.** Fine for the
  single-chunk demos; a world of many chunks needs an addressing scheme that
  names a chunk too. Deferred, not solved.
* **The sink zeroes its own cell, membrane included.** Assembled surfactant *in
  the port cell itself* is not treated as structure. The protection that matters
  happens on the faces between a compartment and the port, which is where it
  belongs, but the port cell is a cruder object than the rest of the model.
* **Outflow is deterministic and exact, therefore mean-field.** A real port
  removes molecules by chance, and stochastic loss is what drives small
  populations extinct — which matters enormously once there are individuals to go
  extinct. Making it stochastic would introduce the engine's first *random*
  physics, a decision Phase 5b deliberately deferred and this phase does not
  pre-empt. Stated, not smuggled.
* **No chemostat demo command yet.** The behaviour is covered by tests; the
  terminal demo is not written.

### New surface

* `au-chem/src/reservoir.rs` — `AtomLedger`, `Port`, `port_delta`, `formula_table`
* `World::atoms` — hashed, snapshotted, restored
* Config: `chem.reservoir.source` (`"cell:species@pop, ..."`), `chem.reservoir.sink` (`"cell; ..."`)
* Tests: +6 (`au-chem` unit) and +11 (`au-sim/tests/openworld.rs`), +2 in `membrane.rs`

---

## Phases 5e, 5f and 5g — changelog reconstructed from source

These three phases shipped as delta archives with their own `PHASE5*_README.txt`
notes and **never had their entries folded into this file**, which is why it
opened at 5c and claimed 164 tests while the tree carried 196. The summaries
below were reconstructed from the module documentation in the source, which is
authoritative for *what the code does*; where the original notes recorded design
history that the code does not, that history is lost and is not invented here.

**Phase 5e — rate constants become physical.** Every generated reaction had
carried the same declared pre-exponential factor, applied to unimolecular and
bimolecular reactions alike — which cannot be right, because the two have
different units (s⁻¹ against m³·molecule⁻¹·s⁻¹). The constant was a fitting
parameter that worked only because a second error cancelled it: the engine ran at
concentrations far below anything real and a prefactor far above anything real
brought the product back into range. Two compensating errors are survivable until
something has to agree with one of them, and membranes arriving with genuine
numbers — a millimolar CMC, a 0.4 nm² headgroup — ended that. `au-chem/src/kinetics.rs`
now derives prefactors from collision theory: Eyring's `k_B·T/h` for unimolecular
steps, `P·σ·√(8k_BT/πμ)` for bimolecular ones with σ and μ read off the
discovered graphs, and an encounter-volume penalty for the termolecular catalysed
reactions of 5c — which is *why* three-body reactions are rare, a fact 5c had
found empirically before there was any physical reason for it in the code.
Forward and reverse now carry different prefactors, and their ratio is the
reaction entropy: a term listed as a known omission for three phases, not added
but falling out of counting molecules on each side.

**Phase 5f — amphiphilicity and self-assembly.** A RAF is a property of a whole
pot: one soup, so nothing is an individual, and with no individuals there is
nothing for selection to act on. `au-chem/src/amphiphile.rs` asks a narrow
structural question of a discovered graph — can it be cut into a polar piece and
an apolar piece? — scoring every **bridge** bond by the same electronegativity
table that prices every bond in the engine, so rings correctly score zero because
a ring has no ends. `au-chem/src/assembly.rs` turns that score into behaviour via
the critical micelle concentration: free monomer rises to the CMC and then stops,
every further molecule joining an aggregate instead, so a structure grows by
accretion rather than dilution — bags, not brine. Nothing is stored: the
assembled material in a cell is a function of the populations already there,
recomputed when needed, exactly as temperature is computed from energy and never
stored. A world's membranes survive save and resume because they were never
separately saved.

**Phase 5g — membranes become selectively permeable.** `diffuse_explicit_perm`
and `diffuse_implicit_perm` add a per-cell permeability that multiplies the
exchange across each face, derived each tick from assembled coverage. A membrane
does not close a face; it slows it. Assembly is also treated as a sink during
transport — a molecule that is part of a structure has no independent gradient to
run down — with the same integers subtracted and restored, so conservation and
non-negativity are untouched. Without that, a compartment would dissolve itself
by diffusion the moment it formed.

## Phase 5c — Catalysis and autocatalytic sets: the chemistry makes itself (incremental update)

**164 tests passing** (was 143). Every earlier test passes unchanged; worlds without `chem.catalysis` (the default) run byte-identically to Phase 5b. This is the third and last prerequisite of Phase 5: 5a let the world invent molecules, 5b let them travel, and 5c is about whether they can make *each other*.

---

### The finding that motivated the phase

Before writing any new physics, the obvious question: how far is the discovered network from self-replication? The answer turned out to be *infinitely* far, for a structural reason.

Form, Split, Raise and Lower give stoichiometries of 2→1, 1→2 and 1→1, and **none of them can place a species on both sides of the arrow**. A molecule that cannot appear among its own reaction's reactants can never make more of itself. Self-replication was not merely absent from the Phase 5a network — it was impossible in it, at any temperature, for any length of time, on any planet. `the_derived_network_cannot_be_stoichiometrically_autocatalytic` proves this over a real 232-reaction network rather than asserting it, and `au raf` reports the same zero in Act I.

Catalysis is what builds the road. A catalyst enters a reaction and leaves unchanged — it appears on *both* sides — and when the catalyst is also the product, the stoichiometry reads `A + B + AB → 2 AB`.

### Catalysis, derived from structure

`chem.catalysis` (default off) enables a **template model**, and it is named that because that is what it is. C catalyses the joining of A and B if it can grip both at once: two free valence slots — which may be two slots on one bridging atom, as the oxygen in A–O–B does — with each grip scored by the same bond-energy model that prices every real bond in the engine. Only *polar* grips count, which is the crude form of the true statement that dissimilar atoms attract.

What it captures: proximity, orientation, and the fact that catalytic ability is a property of structure. What it does not: electronic stabilisation of a specific transition state, geometry, strain, solvent, or anything resembling an active site. Isomerisation (Raise/Lower) has no joined pair and so is never catalysed — a documented gap. The model's virtue is that no molecule was ever declared a catalyst: a species invented on tick 900 is scored by the identical rule as one present at boot.

### The law a catalyst may not break

**A catalyst changes rates. It never changes equilibria.** A catalyst that shifted an equilibrium would let you run a cycle forward over it and back without it and extract work forever. This is the anchor the phase is validated against, and two implementation details are load-bearing:

- The barrier reduction is an **absolute** quantity in J·mol⁻¹, never a fraction of each barrier. Phase 5a built `Ea_f = B + max(ΔH,0)` and `Ea_r = B + max(−ΔH,0)`, so a fixed *fraction* would scale both and change their difference — quietly moving every catalysed equilibrium. That is the kind of error that yields a plausible-looking simulation which is silently a free-energy machine.
- The reduction is capped strictly below the intrinsic barrier `B`, the smaller of the two, so neither direction can be driven to zero and `Ea_f − Ea_r = ΔH` survives every clamp.

### Three errors, all caught, all instructive

**The exclusion rule was asymmetric.** The first `catalysts_for` excluded any species appearing among a reaction's reactants or products. For `A + B → AB` the molecule AB is a product; for the inverse `AB → A + B` it is a reactant — so that rule let AB catalyse the forward direction and forbade it from catalysing the reverse. A catalyst present in one direction only is precisely an equilibrium-shifting machine. The exclusion is now stated over the *joined pair*, which is symmetric under inversion by construction.

**Bridging catalysts were excluded.** The rule demanded two free sites on *distinct atoms*, which rejects the commonest bridge in chemistry — one atom holding two partners. In a world of small molecules that removed nearly every catalyst there was, and the integration tests reported it bluntly as "enabling catalysis changed nothing". Slots are now counted with multiplicity.

**And then catalysis still changed nothing** — which was not a bug at all. At the original concentrations *every* reaction already saturated the reactor's consumption cap, so a faster road was invisible; and a catalysed reaction is termolecular, paying a concentration factor for having to find its catalyst as well as both reactants. A 6× barrier gain cannot cover that. **Catalysis only wins where the barrier it removes is worth more than the encounter costs** — which is exactly why enzymes buy factors of 10⁶ and up rather than 6. The tests and the demo now run in a regime where Arrhenius still gates the rate, with the reasoning written where the config sets it.

### RAF detection: an instrument, not a mechanism

`au_chem::raf` implements Hordijk & Steel's (2004) closure-and-pruning algorithm for the maximal RAF — Reflexively Autocatalytic and F-generated: every reaction catalysed by something the set itself can build from food, every substrate either food or built by the set. Polynomial, deterministic, sets and integers in sorted order.

It changes nothing. `au_sim::autocatalysis::survey` enumerates on a **clone** of the registry so that even looking cannot intern a species, and a test asserts the world hash, clock and registry are untouched. Delete the module and every world runs bit-identically.

Notably, **"a RAF appeared" is deliberately not an event.** An event is world state — hashed, snapshotted, resumed — and the emit-once rule would depend on knowing that the previous tick had no RAF, which no snapshot carries; a resumed world would re-announce a discovery it had already made and diverge from one that never stopped. When a timeline system arrives it can record this from serialized state. Until then it stays a lens.

The **food set** is taken to be the monatomic species — the geochemically meaningful reading, since volcanism and weathering supply atoms. The idealisation is stated: a sealed demo cell does not actually replenish them, so a RAF found there answers "could this network sustain itself if fed?" rather than "is it sustaining itself now". Sustained-in-fact needs flux boundary conditions the chemistry does not yet have.

### What `au raf` shows

Three elements, one hot cell, nothing else declared. Without catalysis: 37 species, 100 reactions, **0** with a species on both sides, no core. With catalysis: 27 catalyst species, a **45-reaction autocatalytic core** closing on 21 species from bare-element food in 3 pruning rounds, and five species that make more of themselves — CO, HCO, C₂O, and two isomers of CO₂ the world discovered separately. The reaction, as the reactor runs it:

```
C + O + CO  →  CO + CO
```

Carbon monoxide holding the two atoms together while a second carbon monoxide walks out. Atom audit exact; save/resume identical.

And **Kauffman's phase transition, in a network nobody designed.** Sweeping how polar a grip must be to count, over one fixed network: 45 reactions with 13 self-producing species → 2 with 0 → nothing. Self-sustaining chemistry does not fade in; it switches on. The cliff sits exactly where carbon–oxygen contacts stop counting — take away C–O polarity and what survives is a vestige with nothing in it that makes more of itself. Nobody put carbon there.

### Honest limits

- **A RAF is not life.** No membrane, no heredity, no individuality; it cannot evolve. It is metabolism before there is anything for the metabolism to belong to — the network-level precondition, not the thing itself.
- The template model is promiscuous: in a rich soup a large fraction of species catalyse something. Specificity beyond polarity is future work.
- `chem.catalysis.max_variants` (default 4096) caps the catalysed set, keeping the strongest reductions. A declared edge like the species pool: deterministic, stated, and reported rather than hit silently.
- Equilibria remain enthalpy-only (no ΔS), inherited from 5a.

### New surface

- `au-chem/src/catalysis.rs` — `CatalysisRules`, `catalytic_reduction`, `catalysts_for`, `catalysed_variants`
- `au-chem/src/raf.rs` — `closure`, `maximal_raf`, `self_producing`, `RafReport`
- `au-sim/src/autocatalysis.rs` — `survey`, `survey_with`, `Survey`
- Config: `chem.catalysis`, `chem.catalysis.min_grip_j_mol`, `.max_reduction_frac`, `.grip_half_j_mol`, `.max_variants`
- `au raf` — the three-act demo
- Tests: +14 (`au-chem/tests/raf.rs`) and +7 (`au-sim/tests/autocatalysis.rs`)


## Phase 5b — Species transport: the chemistry travels (incremental update)

**143 tests passing** (was 128). Every earlier test passes unchanged; worlds with `chem.diffusion_m2_s` unset (the default, 0) run byte-identically to Phase 5a. This is the second prerequisite of Phase 5 (Origin of Life): 5a let the world invent molecules; 5b lets the inventions *go somewhere*. Life will not emerge in a cell that keeps everything it makes.

---

### What changed, in one sentence

Species populations now diffuse between cells by Fick's law — conduction's exact twin, with concentration in temperature's chair — at per-species rates **read off the discovered molecular graphs**, wherever the host material is fluid.

The headline, runnable as `au vent`: hydrogen and oxygen atoms injected at one cell of a molten tube. The world discovers water there — and then water is found at the tube's ends, having been carried by nothing but its own gradient. In deep time (10⁴-second strides, k ≈ 200–700, far past the explicit ceiling of ½) the ocean levels to a **worst deviation of 1 count on a mean of 59,723**, the atom audit stays exact to the last of 6,000,000 atoms, and a world saved mid-spread resumes bit-identically.

### Graham's law: the scale is declared, the shape is derived

`chem.diffusion_m2_s` is the coefficient of a 1-dalton reference particle. Every species then derives its own rate as D·√(1 Da/m) — Graham's 1846 law — with the mass computed from the graph the world drew for itself. A molecule invented yesterday already knows how fast it moves today. To be precise about what is validated: the √m scaling is *input*, not output; what the tests prove is that the tower delivers it unmangled — mass read correctly off the graph, per-species rates threaded through the scheduler, transport exact — with H· outrunning O· by the analytic 3.98 end-to-end. (Ra_c played this role for the fluid solver.)

The `vent` demo measured an H₂/H₂O spread ratio *above* Graham and the reason is visible in H₂'s bimodal profile: **the vent eats hydrogen where the oxidant is concentrated, so the H₂ that survives is the H₂ that escaped.** Reaction–diffusion coupling, unasked for, exactly the kind of thing this project exists to produce.

### Solutes, and the freezing-ocean lever

Species are modelled as solutes in the host material: they move only where the host's phase — derived by the same `derive()` everything else trusts — is liquid or gas. Solids and vacuum are walls. One rule, one consequence worth the whole abstraction: **a cooling ocean that crosses its freezing point traps its chemistry in place**, a frozen stratigraphy released by any later thaw. Tested: a 300 K tube holds its injected species to the exact count, in the exact cell.

### Two bugs, both real, both documented

**The donor's-pocket limiter.** The first one-sweep solve at k = 300 drove populations *negative* — the sum conserved (symmetric transfers cannot lose matter), but impossible physics slipped through. The guard is the finite-volume classic: no face may export more than the donor's holdings divided by its open-face count. Order-independent, exactly conservative, inactive whenever the solve is honest. It makes wrong physics safely wrong instead of impossibly wrong. Its cost is also stated: at extreme k it rate-limits a stride to about one cell of advance per face — stability, not teleportation.

**The deep-time checkerboard.** At k ≈ 1000 the tube's end cells starved to exactly zero while their neighbours held double shares. A seven-cell probe found the mechanism: at large k the implicit system approaches a Laplace problem, where single-grid Gauss–Seidel's convergence factor is set by the *grid*, not by k — about cos²(π/n) per sweep — so 24 sweeps left rows violated by ~10⁶ counts; the transport multiplied that error by k, manufactured a red-black checkerboard, and the warm start re-imprinted it every stride. The fix is not cleverness but honesty about the method: **sweeps must scale as 2·n²** (`sweeps_for_deep_time`), after which the same 41-tube mixes to the *exact* uniform count with residual 0.0 — locked by a regression test. Conduction's solver shares this trait; its validated regimes sit below the threshold, its transport divides by heat capacity (blunting the amplification), and unifying both under multigrid is noted as the future fix.

### Honest limits

- **The quantisation deadband is the spec.** Transfers round to integers, so gradients shallower than ~1/(2k) counts per face persist; the equilibration test pins "no two neighbours differ by more than 1" as the terminal state. The honest microphysics below that floor is stochastic Brownian hopping, which would make diffusion the engine's first *random* physics — a real decision about the determinism story, deferred and documented rather than smuggled in via rounding tricks.
- **Advection is not here.** Species ride their gradients, not the currents; coupling them to the Boussinesq velocity field is a natural follow-on.
- **Chunks remain sealed to each other**, for species exactly as for heat — one halo mechanism will serve both, later.
- **The explicit path substeps** to the stability ceiling and refuses pathology loudly: past 10,000 substeps per tick it panics with instructions to set `physics.implicit_conduction = true` — the same one flag that takes heat and matter into deep time together.

### New surface

- `au-physics/src/diffuse.rs` — `diffuse_explicit`, `diffuse_implicit`, `DiffuseScratch`, `sweeps_for_deep_time`
- `chem.diffusion_m2_s` — the one new dial (default 0 = off)
- `au vent` — the two-act demo
- Tests: +10 (`au-physics/tests/diffuse.rs`: exact conservation at any convergence, Einstein's σ² = 2kt within 3%, implicit/explicit agreement, walls, the deep-time regression) and +5 (`au-sim/tests/diffusion.rs`: Graham end-to-end, frozen trapping, deep time, discovered water spreading beyond its birthplace, save/resume mid-spread)


## Phase 5a — Open-ended chemistry: the engine discovers (incremental update)

**128 tests passing** (was 115). Every earlier test still passes unchanged; closed, declared chemistry runs byte-identically. This is the first prerequisite of Phase 5 (Origin of Life), and it removes a ceiling that has been in place since Phase 3: **a world whose molecule list is written in config can never contain a surprise.** Life cannot emerge from a closed vocabulary. Now the vocabulary is open.

---

### What changed, in one sentence

Reactions are no longer declared; they are **derived from molecular structure** by local graph moves, and when a move produces a molecule the world has never seen, it is interned on the spot — discovered, logged as history, folded into the world hash, serialized in every save.

The headline, runnable as `au genesis`: a config that declares **two kinds of atom and nothing else** — hydrogen and oxygen — is put in a hot sealed cell. Within five ticks the world has discovered H₂, hydroxyl, O₂, water, and HO₂. Nobody listed water. It appears because it is reachable, and it wins because it is the deepest well in the neighbourhood.

### The moves, and why they come in inverse pairs

A candidate reaction is one application of one move: **Form** (a new single bond between free-valence atoms of two molecules), **Split** (remove a single bond whose loss disconnects), **Raise** (bond order +1 where both atoms can pay), **Lower** (the opposite). Form/Split are exact inverses; so are Raise/Lower. That pairing is thermodynamics, not tidiness — every transformation the generator can produce, it can also reverse, so the network settles into genuine equilibria instead of one-way ratchets. Deliberately excluded, each documented at the code site: ring closure/opening (they must arrive as a pair, and the H/O validation world does not need them) and concerted substitutions (representable as split-then-form; collapsing them changes kinetics, not reachability).

### The activation identity everything rests on

Generated reactions need activation energies, and the model is the simplest one that preserves the thermodynamics Phase 3 validated: `Ea_f = B + max(ΔH, 0)`, `Ea_r = B + max(−ΔH, 0)`, one declared intrinsic barrier `B` for the whole chemistry. Subtracting gives `Ea_f − Ea_r = ΔH` **exactly, by construction** — so with equal pre-exponentials, `K = exp(−ΔH/RT)`: Boltzmann, structural. A test pins the identity across the entire generated closure, because "true by construction" is a property of code and code gets refactored. A second test runs a reactor on purely *generated* reactions (2H· ⇌ H₂, both directions found by the machine) and measures the equilibrium at two temperatures against the Boltzmann prediction: within 25% at 3000 K and 4500 K. Chemistry the engine invented obeys laws nobody wrote into the generator.

### Why enumeration does not explode

Attaching hydrogen to any of benzene's six carbons gives one product, and enumerating it six times would be waste. The generator works on **symmetry classes** — the converged Morgan-refinement ranks, now exposed as `Molecule::symmetry_classes()` (extracted from `canonical()`, which is unchanged and untested-code-identical). One representative site per class for Form; one representative bond per class for Split/Raise/Lower. Water offers one O–H bond to split, not two. Because refinement can in principle over-merge on regular graphs, finished candidates are *also* deduplicated by the canonical forms of their full term lists — the classes keep enumeration small, the canonical dedup keeps it correct, and neither is trusted alone.

### Shelf space, not a cliff

`chem.max_atoms` bounds molecule size (the combinatorial guard), and `chem.max_species` bounds the world's column pool. Discovery past the pool is *recorded but never stocked*: the molecule enters the registry — the world knows it — but reactions touching it are dropped from the active set, so no population can ever exist outside the columns that were registered at boot. A declared edge, stated rather than hit; a test crams a world into a pool of three and demands determinism rather than a crash.

### The registry became history — and the save format grew a frame

Under declared chemistry the species registry was a pure function of config, rebuilt at boot, never saved. Discovery changes its nature: which molecules exist is now something this world *did*. So the registry serializes (`SpeciesRegistry::to_bytes`, molecule graphs in id order), rides in every save, and folds into the world hash — two worlds that found different molecules are different worlds before a single population differs. Each interning also emits a `CHEMISTRY_DISCOVERY` event (id, atom count, total), because the moment a substance began to exist is a cause like any other.

The save format is now a frame: `[u64 snapshot length][snapshot][registry]`. The au-data codec is untouched — the snapshot inside is byte-for-byte the format it was — but **saves from earlier phases no longer load**, and this line is the notice. No compatibility machinery was added because no compatibility promise exists yet; when one does, the frame gives it a place to live.

### What the demo surfaced, unasked: kinetic trapping

The genesis run produced a finding the code was never asked for. At 1500 K, water is discovered in the first five ticks — and then almost none is made. The fast reactions run first: every free O grabs an H (hydroxyl, ~2,000,000), leftover H pairs off (H₂, ~1,000,000), and the cell goes quiet with 58 molecules of water, because the road onward runs through cracking H₂ — a 250 kJ/mol toll that 1500 K pays almost never. The matter sits in the shallow well it reached first, not the deep one it can see. Heat the same soup to 3000 K and the barriers open: after 600 ticks, **99.86% of the oxygen is water.** Kinetic trapping and thermal annealing, emerging from four graph moves and one activation identity — real chemistry's signature behaviour, and precisely the kind of unprogrammed consequence this project exists to produce. (It is also why metastable things — proteins, planets, readers — exist at all.)

### Honest limits, stated

One cell: transport of discovered species between cells is Phase 5b, and its absence here is a scope line, not an oversight. Per-event reaction heat is now at the physically honest scale (molar enthalpy / N_A ≈ 10⁻¹⁹ J), so genesis barely warms its cell — Phase 3b already proved the heat coupling; genesis is about matter. The equilibrium model is enthalpy-only (`K = exp(−ΔH/RT)`, no ΔS term), which over-favours association at low temperature relative to real free energies; an entropy model is a refinement with a place reserved. `formula_string` is cosmetic output, read by nothing.

### Files

New: `crates/au-chem/src/network.rs`, `crates/au-chem/tests/network.rs` (9 tests), `crates/au-sim/tests/genesis.rs` (4 tests). Changed: `crates/au-chem/src/{lib.rs, molecule.rs, reaction.rs}` (exports; `symmetry_classes` extraction; registry serialization), `crates/au-core/src/event.rs` (`CHEMISTRY_DISCOVERY`), `crates/au-sim/src/{chemistry.rs, world.rs, systems/chemistry.rs, sim.rs}` (open config; registry on `World`, seeded, hashed; the open-mode system; the framed save/load), `crates/au-headless/src/main.rs` (the `genesis` demo), `README.md`, `CHANGELOG.md`.

---

## Phase 4b — Implicit conduction (incremental update)

**115 tests passing** (was 103). Every earlier test still passes **unchanged**, and every earlier phase reproduces bit-for-bit: implicit conduction is off by default. This phase closes the project's oldest structural crack.

---

### The crack

Every phase since Phase 2 has been honest about a limit it could not pass. An explicit diffusion scheme is stable only while `dt ≤ dx²/(2·ndim·α)`, and that `dx²` is fatal — halve the cell size for detail and the timestep must quarter, so a fixed span of simulated time costs the *fourth* power of resolution. Phase 1's deep-time accelerator turned out to be largely inert against real physics for exactly this reason. Weather needs hours, geology needs epochs, evolution needs generations, and none of them were reachable. Simulating one million years of rock at metre resolution needs ~3×10⁸ explicit substeps — for a single chunk.

Rather than build Phase 5 upward past a known crack, this phase goes down into the Physics layer and fixes it.

### The fix

Backward Euler. Instead of computing the new field from the old one, define it in terms of itself — `(I − dt·α·∇²)·T^{n+1} = T^n` — which turns one explicit sweep into a sparse linear system that is diagonally dominant and nearest-neighbour: the same shape as the pressure projection Phase 2b already solves, so the red-black SOR built there solves this too. The reward is that **no timestep is unstable**. A step a thousand times past the explicit limit is merely less accurate, which is a trade a caller can make knowingly.

Measured: 5000× fewer substeps and ~800× less wall clock for the same span of physics, with the advantage growing without bound as the grid is refined.

### An approximate solver inside an exact engine

The real difficulty is not the solve; it is that an iterative floating-point solver converges to a *tolerance* while the whole project rests on energy being an exact `i128`. Bolted in naively this puts a rounding leak at the centre of the foundation — and a leak is exactly what a later evolutionary optimiser learns to farm.

The resolution is to give the two pieces different jobs. The **solver** decides *how much energy should move*: approximate, iterative, allowed to be imperfect. The **transport** moves it: one integer per face, subtracted from one side and added to the other. So **conservation does not depend on the solve converging.** An under-converged solve puts heat in the wrong place — wrong physics, honestly wrong, and visible in the reported residual — but it cannot lose a joule, because losing one would require the two halves of a symmetric integer transfer to disagree, and they are the same integer. A test runs the solver deliberately crippled at one sweep and demands the books still balance to zero.

This is the third time the project has needed the same idea: Phase 2's symmetric flux, Phase 3's "quantise the sum, never the summand", and now this. Approximation may make the answer wrong; it may never make it non-conservative. The corollary is a trap worth naming, because it is the obvious implementation: do **not** write the update as `ΔE = C·(T^{n+1} − T^n)` per cell. That is conservative only if the solve converged perfectly, silently coupling the central invariant to an iteration count.

### Two real bugs, both caught by the tests

- **A sentinel that grew a second meaning.** `stable_dt` returns infinity when nothing constrains the timestep. That used to have exactly one cause — an inert grid with no diffusivity — so an early return was a harmless shortcut. Implicit conduction made infinity *also* mean "the stiff term is unconditionally stable, take the whole step", and the old shortcut therefore skipped the entire physics step. Every implicit run was a silent no-op. Tellingly, the conservation tests passed happily throughout — doing nothing conserves beautifully — and only the physics tests (equilibration, downhill flow) caught it. Fixed by taking one full stride, which is correct for both meanings.

- **Over-relaxation forfeits the maximum principle.** The first version passed the pressure solver's `sor_omega` of 1.85 to the heat solve, and the stability test found temperatures of 4244 K in a box that started between 400 K and 1200 K — energy still conserved to the microjoule, but heat piled up where physics forbids. The cause is structural and worth writing down: the conduction update is a weighted average of a cell's own previous temperature and its neighbours', with positive weights summing to one — a *convex combination*, which cannot leave the range it started from however few sweeps run. Over-relaxation breaks the convex combination, and because the transported energy is that temperature difference times a very large `dt`, any overshoot is amplified rather than damped. `conduction_omega` therefore defaults to 1.0 (plain Gauss–Seidel): slower to converge, but the discrete maximum principle holds even on an absurdly under-converged solve. A test pins that at one sweep and a step a million times past the explicit limit.

### What stays explicit, and why

Radiation is nonlinear in `T⁴` and would need Newton iteration; its bound scales as `T³`, not `1/dx²`, so it does not worsen under refinement. Advection is bounded by `dx/u`, linear in `dx`. Viscosity has the same `dx²` problem as heat and deserves the same cure — it is the next candidate, and this module is written to be reusable for it. Stability is now unconditional; **accuracy is not** — a huge step still smears fast transients, because backward Euler damps what it cannot resolve. The cliff is gone; the hills remain.

### Files

New: `crates/au-physics/src/implicit.rs`, `crates/au-physics/tests/implicit.rs` (12 tests). Changed: `crates/au-physics/src/lib.rs` (module + exports), `crates/au-physics/src/grid.rs` (`implicit_conduction`, `conduction_iters`, `conduction_omega` on `Boundary`), `crates/au-physics/src/transport.rs` (scratch field, the `stable_dt` exemption, the conduction route, and the sentinel fix), `crates/au-sim/src/systems/physics.rs` (config keys `physics.implicit_conduction` / `.conduction_iters` / `.conduction_omega`), `crates/au-headless/src/main.rs` (the `deeptime` demo), `README.md`, `CHANGELOG.md`.

---

## Phase 4 — Planet generation (incremental update)

**103 tests passing** (was 84). Every earlier test still passes **unchanged**. New crate `au-planet` (14 tests), a planetary chunk generator in `au-sim` (5 integration tests), and the `au planet` demo. From this phase on, a world can declare a planet in config and its chunks fill themselves — the first matter in the project that exists without a test fixture placing it.

---

### The fork, and the road taken

There are three ways to make a planet. **Noise terrain** (a Perlin heightmap painted with biomes) violates the Golden Rule outright — the mountains would exist because a noise function put them there, a river drawn rather than carved. **Full physical formation** (accretion to tectonics, dust to world) is the honest dream and an uncomputable cliff — and Phase 2 already proved deep time cannot be fast-forwarded through diffusion, so we could not reach 4.5 billion years even in principle. The road taken is the third: **a planet is *declared* by its physical parameters** — mass, radius, orbital distance, stellar luminosity, albedo, rotation, surface pressure — **and its gross state is *derived* by the physics already built.** Gravity from Newton. Equilibrium temperature from Stefan–Boltzmann — the very law Phase 2 validated against Stefan (1879), now used globally. Ocean-or-ice-or-vapour from water's own `melt_at`/`boil_at` phase model — the identical code that melts an ice cube in a Phase 2 cell. Atmospheric scale height from hydrostatics. Phase 4 adds almost no new physics; it aims the existing physics at a planet.

This is the exact analogue of the chemistry decision: declare the vocabulary (which elements exist / how massive, how far), derive the outcome (which reactions fire / whether water stands liquid). Nobody paints an ocean. The star's warmth and water's boiling point decide there is one.

### The validation: a band nobody drew

Phase 2b had Ra_c = 657.5; Phase 3 had van 't Hoff. Phase 4's prediction is the **habitable zone**. Sweeping orbital distance, the surface volatile goes vapour → liquid → frozen, monotonically outward, in a single contiguous band — and a 4× brighter star pushes that band outward by the √L the inverse-square law demands. Earth's parameters yield 9.82 m/s², 1361 W/m², and 255.5 K — the textbook effective temperature — plus a 7.5 km atmospheric scale height. All derived, none stored.

Two honest limits surfaced by the physics itself, both now documented where they bit:

- **The bare-rock Earth freezes.** 255 K is below water's freezing point; the missing ~33 K is the greenhouse effect, and this phase's atmosphere is transparent to its own thermal radiation. Even a perfectly black Earth (albedo 0) only reaches ~277 K — you cannot close the greenhouse gap with albedo, because the gap *is* the greenhouse. Tests use an albedo-0 world as the declared stand-in, and say so.
- **Venus is hell because of its atmosphere, not its orbit.** At 0.72 AU the bare radiative balance is only ~345 K — below boiling. Without a runaway greenhouse, a planet at Venus's distance would hold a hot but liquid ocean. My first vaporisation test assumed otherwise and the engine correctly failed it; the test moved to 0.40 AU, where radiation alone exceeds boiling, and the lesson moved into a comment.

### Terrain: the licence, written down

Relief cannot yet be caused (orogeny is diffusion's slowest cousin — it needs the implicit viscous solver, the same one circulation needs). So relief is **initial-condition detail**, the same epistemic category as the seed: a declared amplitude and wavelength, the particular highs and lows drawn from position-keyed value noise (`Rng::derive(seed, PLANET, x, y)` — the splittable-RNG discipline, which also makes the surface seamless across chunks generated independently). The noise decides only *where the rock stands a little higher*. Everything about what that rock **is** — temperature, phase, whether the hollow beside it floods — is physics. When tectonics arrives, relief becomes an output of simulation and the noise retires to seeding the primordial surface. The boundary is documented in `terrain.rs`, at the spot where the dice are rolled.

### The generator: seed decides where, physics decides what

`PlanetGenerator` (in `au-sim`, where planets meet columns, as chemistry did) fills each cell by two questions with two different authorities. *Where does the rock stand?* — the terrain field answers. *What is everything's state?* — physics answers, and only physics. The generator writes mass and energy; the energy implies the star's derived temperature; the **phase is whatever `derive` reports when anyone looks**. The headline integration test runs the same config at two orbits: at 0.95 AU the below-sea-level cells derive as **liquid** — an ocean; at 1.6 AU the *same cells from the same code* derive as **solid** — an ice sheet. The generator never chose. Five integration tests: the world fills its own chunks (no hands), columns are causally ordered (rock, then volatile, then vacuum), warm-makes-ocean/cold-makes-ice, bit-determinism, and save/resume — where `Simulation::load` now reinstalls the planet generator automatically, because a planet is config and config rides in the snapshot.

### The demo

`au planet` derives Earth's numbers live, renders the habitable zone as a band (ocean 0.56–1.04 AU for the albedo-0 proxy), then generates the same seed at 0.95 AU and 1.6 AU side by side: identical relief, sea in one cross-section, ice sheet in the other, both worlds reloading hash-identical.

### Deliberately deferred

Atmospheric circulation, weather, and hydrology (all need the implicit solver); greenhouse radiative transfer; accretion, differentiation, tectonics; per-cell volumes and partial-cell coastlines; eccentric orbits and seasons; deriving surface pressure from outgassing and escape. Each deferral is written at the code site it constrains.

### Files

New: `crates/au-planet/` (Cargo.toml, src/{lib,planet,terrain}.rs, tests/validate.rs), `crates/au-sim/src/planet_gen.rs`, `crates/au-sim/tests/planet.rs`. Changed: workspace `Cargo.toml` (member), `crates/au-core/src/rng.rs` (a `PLANET` RNG domain), `crates/au-sim/Cargo.toml` + `src/lib.rs` + `src/sim.rs` (dependency, module, auto-install in both boot and resume), `crates/au-headless/Cargo.toml` + `src/main.rs` (the `planet` command), `README.md`, `CHANGELOG.md`.

---

## Phase 3b — Chemistry in the simulation loop (incremental update)

**84 tests passing** (was 79). Every earlier test still passes **unchanged**. This phase wires the Phase 3 chemistry engine into the real simulation — it runs in the scheduler, on the active set, feeding the same thermal field physics conducts and radiates — without touching the engine's guarantees.

---

### What Phase 3 left open, and this closes

Phase 3 built and validated the chemistry engine in isolation: 18 tests, van 't Hoff reproduced, exothermic reactions heating a cell with exact conservation. But it sat *beside* the simulation, not in it. Phase 3b puts it in the loop:

- A **chemistry system** (`systems/chemistry.rs`) runs each tick after physics, over the active set, following the same gather → compute → scatter shape the physics system uses.
- It reads and writes the **same `PHYS_ENERGY` column** the thermodynamics solver owns. An exothermic reaction's heat lands in the field physics then carries away; the resulting temperature feeds back into the Arrhenius rates. One shared, conserved energy account — which was the whole point of building chemistry *on* physics.
- Temperature for the rates is derived from that shared energy by the **same function physics uses** (`au_physics::derive`), so the two layers never disagree about what temperature an energy means.

### The scope decision (and why it isn't cheating)

A world **declares its chemistry vocabulary at boot** — its elements, the molecules that *can* exist, and the reaction rules that *may* fire — in config, exactly the way the material table and the (coming) periodic table are declared. The species count is therefore fixed at boot, which lets per-cell populations live in fixed columns (`0x2xxx` range) that round-trip through the existing snapshot codec with **no format change**.

This declares the *vocabulary*, never the *outcome*. Which reactions fire, how far, how hot a cell gets, whether it ignites — all emergent from rules interacting, every run. Declaring that carbon and water may exist is the same kind of act as declaring which elements the periodic table holds; it is the opposite of scripting a wolf. Open-ended molecular discovery (a reaction inventing a graph no config named) is real and coming, but it needs the species registry to become serialized state and per-cell storage to become sparse — correctly a Phase 5+ problem. The layers above won't change when it arrives; they see species ids and populations either way.

### What the integration tests prove

Beyond "reactions happen" (Phase 3 showed that), the five new `au-sim` tests prove the reactions happen *in the world*:

- **Chemistry heats the shared thermal field in the sim loop** — with boundary flux and radiation off, the only thing that can raise the cell's energy is a reaction, and it does.
- **Atoms are conserved through the coupled run** — H and O only rearrange between H₂/O₂/H–O.
- **Determinism holds** — same seed, bit-identical world hash, now including chemistry's events and every dirty chunk's species populations.
- **Save/resume is identical** — a world saved mid-reaction, reloaded, and finished matches one that never stopped. This is the test that proves chemistry state round-trips through the snapshot codec.
- **A dormant chemistry is inert** — a world with chemistry declared but no reactants present changes nothing and adds nothing to history, so worlds that never react are unaffected by the layer existing.

### The demo

`au react-sim` boots a full `Simulation` with a declared H₂ + O₂ chemistry, activates a cell, and lets the scheduler run physics and chemistry together. The shared energy field climbs as bonds become heat, atoms balance, and the run ends by saving, resuming, and confirming the world hash is **identical** — determinism and resume demonstrated live, not just asserted in a test.

### Files

New: `crates/au-sim/src/chem_columns.rs`, `crates/au-sim/src/chemistry.rs`, `crates/au-sim/src/systems/chemistry.rs`, `crates/au-sim/tests/chemistry.rs`. Changed: `crates/au-core/src/event.rs` (a `CHEMISTRY_HEAT` event kind), `crates/au-sim/Cargo.toml` (depends on `au-chem`), `crates/au-sim/src/lib.rs` (new modules), `crates/au-sim/src/world.rs` (register chemistry columns at boot), `crates/au-sim/src/sim.rs` (register the chemistry system in both the `new` and `load` schedules), `crates/au-sim/src/systems/mod.rs`, `crates/au-headless/src/main.rs` (the `react-sim` demo), `README.md`.

---

## Phase 3 — Chemistry (incremental update)

**79 tests passing** (was 61). Every earlier test still passes **unchanged**. The new crate `au-chem` is additive; nothing in `au-core`, `au-data`, `au-physics`, or `au-sim` changed except that `au-headless` gained a `react` demo and a dependency on `au-chem`.

---

### The honest framing

Real chemistry is not computable at the scale this project runs at. Solving the Schrödinger equation for molecules costs roughly O(N⁷) and runs on supercomputers for *single* reactions; this simulation will have millions of cells. So this crate does not simulate chemistry. It provides an **abstraction of chemistry chosen to preserve the causal structure that lets chemistry produce life**, while staying computable. Getting that abstraction right is the actual work of Phase 3, and everything above it rests on it.

### The three commitments

1. **Elements are exactly conserved integer counts.** Atoms of a given atomic number are rearranged, never created — because a chemistry that can leak a carbon atom is a chemistry that something will evolve to farm. Counts are `i128`, transferred by the same symmetric-flux discipline that got Phase 2 energy to zero error.

2. **Molecules are graphs, discovered not listed.** Water is not a table entry — it is "two H bonded to one O", a labelled graph the engine builds and fingerprints. A canonical form (Morgan-style refinement with McKay-style orbit-breaking for symmetric graphs) recognises when two cells contain the same molecule regardless of how the atoms were enumerated. Genuine novelty is possible because a new molecule is just a new graph.

3. **Reactions are physics-gated rules sharing one energy field.** A reaction fires by mass-action kinetics at an Arrhenius rate, and its enthalpy — *derived* from bond energies, so it cannot secretly violate conservation — flows into the **same** thermal field Phase 2 conducts and radiates. An exothermic reaction heats the cell it occurs in; the warming raises every reaction's rate; and that feedback is combustion. Chemistry and thermodynamics are one conserved loop, not two systems bolted together.

### The result

The engine reproduces **van 't Hoff** to several significant figures. A reversible reaction, left alone, settles at the equilibrium its temperature dictates — and for an exothermic reaction the equilibrium constant *falls* as temperature rises, shifting the balance back toward reactants (Le Chatelier). Nothing in the reactor knows those laws; it only fires rules at Arrhenius rates. Equilibrium is an emergent fixed point, exactly as in reality.

```
temp K     K = [HO]²/[H₂][O₂]
   500                  0.999
   900                  0.962      ← heating an exothermic reaction makes less product
```

And in the `react` demo, a cold parcel of H₂ + O₂ with no external heat source ignites and burns to completion above a threshold temperature — heating *itself*, entirely from H–O bonds forming — while below the threshold the barrier holds. The heat the reactor reports released matches the thermal field's gain to the microjoule: **conservation error 0 µJ**.

### Bugs caught and fixed (all real, all documented in code)

- **Bond-energy sign inversion.** The first model made stronger bonds *higher* energy, which made forming water endothermic — backwards. A molecule with strong bonds sits *lower* in its potential well; fixed so exothermic means products-below-reactants.
- **Kinetic stiffness.** Fast reactions consumed 100% of their reactant per step and oscillated across equilibrium instead of settling — the forward-Euler instability of a stiff ODE. Fixed with a per-step consumption cap and internal sub-stepping, the chemical echo of Phase 2's CFL condition.
- **The dynamic-range floor, again.** A single molecular reaction releases ~5×10⁻¹⁹ J, far below the microjoule quantum the thermal field is stored in, so per-event enthalpy rounded to zero and combustion released no heat. Fixed properly: the reactor carries per-event enthalpy in full-precision joules and quantises only the *aggregate* `events × enthalpy`, carrying the sub-microjoule remainder so nothing is lost. Quantise the sum, never the summand — the same lesson as the Boussinesq mass-quantum bug in Phase 2b.

### What is deliberately deferred to Phase 3b

The chemistry engine is validated in isolation (18 tests) and demonstrated (`au react`). Wiring it into `au-sim`'s scheduler — per-cell species storage in the `0x2xxx` column range, a chemistry system running on the active set at `Layer::Chemistry`, snapshot-format extension, and the determinism work that comes with new persistent state — is substantial and is deferred so it can be done carefully rather than rushed. The scientific claim of Phase 3 is already proven; what remains is integration.

### Deliberately *not* in Phase 3 (per the abstraction bet)

No 3-D molecular geometry, no stereochemistry, no bond angles, no resonance — connectivity plus bond order only. If Phase 5 shows chirality is load-bearing for replication, geometry gets added to `molecule.rs` without any layer above changing, because they only ever see graphs and canonical forms.

### Files

New: `crates/au-chem/` (Cargo.toml, src/{lib,element,molecule,reaction,reactor}.rs, tests/{foundation,reactions}.rs). Changed: workspace `Cargo.toml` (added member), `crates/au-headless/Cargo.toml` (added dependency), `crates/au-headless/src/main.rs` (added `react` command), `crates/au-physics/src/fluid.rs` (removed a no-op `drop`), `README.md`.

---

## Phase 2b — Fluid Dynamics & Convection (incremental update)

**61 tests passing** (was 52). Every earlier test still passes **unchanged** — Phase 1's `universe.rs`, Phase 2's `thermodynamics.rs`, and the 12 thermodynamics tests in `physics.rs` were not touched. A grid with no velocity field is Phase 2 exactly as it was.

---

### The result

The engine reproduces the onset of Rayleigh–Bénard convection. A fluid layer heated from below sits still and conducts below a critical Rayleigh number, and spontaneously overturns above it — and that threshold, for stress-free boundaries, is a number **Rayleigh derived on paper in 1916**:

```
Ra_c = 27π⁴/4 = 657.511…
```

Nothing in the codebase knows that number. The solver is asked only "did this disturbance grow?", and bisection on the answer finds the threshold. On a 12-cell-deep grid it lands at 600; refine the grid and it walks home:

```
cells deep      Ra_c found       error
        12           600.7       -8.6%
        16           610.6       -7.1%
        20           637.2       -3.1%
   Richardson → ∞:      ~684        (brackets 657.5)
```

That convergence *is* the proof. A solver with a bug converges to the wrong number or to nothing; a correct solver with a coarse mesh converges to the analytic answer as the cells shrink. This one walks toward 657.5.

What appears above the threshold is a **dissipative structure** — order, spontaneously organised, sustained entirely by energy flowing through, and gone the instant the flow stops. Nobody places the convection rolls. They are the fluid's answer to a gradient it cannot conduct away fast enough. That is the same *category of thing* life is, which is the whole bet of the project.

### New files

| File | What |
|---|---|
| `crates/au-physics/src/fluid.rs` | The MAC staggered-grid solver: advection, viscosity, buoyancy, and the pressure projection. |
| `crates/au-physics/tests/fluid.rs` | 9 tests — Rayleigh–Bénard onset, Nusselt transport, momentum/mass/energy conservation, incompressibility, spurious-current checks, determinism. |

### Modified files

| File | Change |
|---|---|
| `crates/au-physics/src/quantity.rs` | `Momentum` (exact `i128`), and an impulse ledger — the books that make "a fluid cannot push against nothing" checkable. |
| `crates/au-physics/src/material.rs` | Viscosity (per fluid phase; solids are rigid, not thick) and thermal expansion α — the single coefficient that makes buoyancy, and therefore convection, exist. |
| `crates/au-physics/src/grid.rs` | `Wall` (periodic / no-slip / free-slip), Dirichlet temperature walls, and the staggered `FluidView` — velocities on faces, not at centres. |
| `crates/au-physics/src/transport.rs` | Energy advection (upwind, flux-form); the advective (Courant) and viscous stability limits, stacked on the thermal one; fixed-temperature walls. |
| `crates/au-physics/src/lib.rs` | Exports. |
| `crates/au-physics/tests/physics.rs` | Updated for the new `Material`/`Boundary`/`GridView` signatures. Physics unchanged. |
| `crates/au-data/src/codec.rs` | `FORMAT_VERSION` 2 → 3. |
| `crates/au-data/src/snapshot.rs` | Conserved-quantity ledger widened to 8 slots (momentum impulse joins energy and mass). |
| `crates/au-data/tests/memory_strategy.rs` | Mechanical: the ledger array grew. |
| `crates/au-sim/src/physics_columns.rs` | Momentum columns (0x1004–6) and dynamic pressure (0x1007). |
| `crates/au-sim/src/config.rs` | `f64_opt` (absent ≠ zero — a missing wall temperature means "no wall"), `str_or`. |
| `crates/au-sim/src/sim.rs` | Snapshot v3 wiring for the impulse ledger. |
| `crates/au-sim/src/systems/physics.rs` | Gather → compute → scatter for seven columns; fluid boundary wiring. |
| `crates/au-headless/src/main.rs` | The `au convect` demo, and `au ra --ny N` for the convergence study. |
| `data/reference_materials.kv` | Magma (basalt) gains viscosity/α; a Prandtl-1 calibration fluid added — an instrument, not a substance. |
| `data/constants.kv` | Fluid config: walls, projection iterations, the Boussinesq note. |

### Two decisions, made in the open

**Boussinesq.** Density is treated as constant everywhere except where it multiplies gravity. That is not a detail of the approximation — it *is* the approximation, and it is the one under which Ra_c = 657.511 is a true statement. Fully compressible flow would tie the timestep to the speed of sound (~2 ms in rock) and four billion years would never arrive.

**MAC staggered grid.** Co-locating velocity and pressure decouples odd and even cells in the pressure equation, producing an invisible checkerboard mode that drives flow indistinguishable, to the eye, from real convection. Putting velocities on cell faces (Harlow & Welch, 1965) makes the pressure Laplacian compact and the checkerboard impossible. It makes every loop fiddlier and it is not optional.

### The bug that nearly fooled me

The first solver produced convection at *half* the critical Rayleigh number — with the perturbation set to **zero**. It stirred itself out of nothing. Two causes:

1. **Boussinesq says ρ is constant, and I was advecting mass anyway.** Under a divergence-free field that is exactly zero in the continuum, so I was integrating pure rounding error and calling it transport.

2. **The signal was below the quantum.** Temperature is derived as E/(m·c). The mass quantum is 1 µg in a 10⁻³ kg cell — a relative precision of 10⁻⁶ — while my driving ΔT was 6.5×10⁻⁷ relative. I had asked the engine to resolve something finer than the grain of its own state, and it answered with noise, correctly.

This is the price of exact integers, and it is worth naming: they buy perfect conservation and charge for it in **dynamic range** — a hard, knowable floor beneath which nothing is real. Floating point would not have complained; it would have leaked energy *and* produced the same fictitious convection, and I would have believed it. The fix: Boussinesq does not advect mass, and the calibration is chosen so criticality sits at ΔT ≈ 1 K, three thousand times the quantum. The test that would have caught it now exists: `a_stratified_subcritical_fluid_does_not_stir_itself`.

### New conservation law

**Σ momentum(now) − Σ momentum(start) == impulse**, exactly. Advection and viscosity only move momentum between control volumes (one integer out of one, into the other); everything external — buoyancy, the pressure force on the walls, wall friction — is booked. A fluid that gains momentum from nothing is a reactionless drive, and to a selection process free propulsion is the same discovery as free food. Demonstrated at `0` in the demo, with the fluid in full flight.

### Performance

The pressure solve sweeps the grid 60× per substep. The first version recomputed cell coordinates and wall tests inside that loop — three-quarters of total runtime. Precomputing the Poisson stencil once per substep and looking it up: **10× faster**, and it is why the demo runs in 15s instead of minutes.

### Not done, on purpose

- **No 3-D convection cells.** The solver is fully 3-D and conserves in 3-D; the demo is a 2-D slice because it is legible in a terminal and because Ra_c is a 2-D result. Nothing structural blocks 3-D.
- **Real mantle viscosity is unreachable by this method.** The mantle's kinematic viscosity is ~10¹⁷ m²/s, which demands a viscous timestep of ~10⁻¹³ s. Explicit integration cannot go there. Reaching geological time needs an implicit viscous solve — a Phase 4 concern, and now a known one.

---

## Phase 2 — Thermodynamics (incremental update)

**52 tests passing** (was 31). Every Phase 1 test still passes **unchanged** — `crates/au-sim/tests/universe.rs` was not touched.

---

### The decision that shaped this phase

Physics for this project is **thermodynamics and transport**, not rigid-body mechanics.

Nothing above this layer needs a falling crate. Chemistry needs temperature, pressure and available energy. Climate *is* heat transport. And the origin of life needs a **persistent free-energy gradient** — life is a dissipative structure; it grows in the gap between a hot thing and a cold thing. Building a collision engine would have been building the wrong thing beautifully.

### The constraint that made it hard

**Evolution is an adversarial optimiser whose search space includes the bugs in your physics.** Leave any path by which energy can be created and something will eventually evolve to walk it — not because it is clever, but because natural selection is an exhaustive search and free energy is the most rewarding thing in any possible fitness landscape. This is the classic way artificial-life projects die.

So conservation is not "accurate to 1e-9". It is **exact**:

* Energy (µJ) and mass (µg) are `i128` integers. They cannot drift.
* Transport is **symmetric flux** — one integer, subtracted from one cell and added to the other. Conservation is an identity, not an aspiration.
* Floats appear only in **derived** values (temperature, density, pressure), recomputed fresh each tick and never accumulated. Float error cannot compound if you never compound floats.

Measured conservation error over 2,000,000 ticks with heat flowing in and out: **0 µJ**.

### New files

| File | What |
|---|---|
| `crates/au-physics/Cargo.toml` | New crate. Depends only on `au-core` — it knows nothing of worlds, chunks or planets, which is what lets it be tested against analytic solutions with no simulation running. |
| `crates/au-physics/src/quantity.rs` | `Energy`, `Mass`, `Ledger`. Exact integers, and the argument for why. |
| `crates/au-physics/src/material.rs` | Bulk substances: heat capacity, conductivity, latent heat, Clausius–Clapeyron slopes. **Registry ships empty** — Phase 3 will *derive* materials from elements rather than list them. |
| `crates/au-physics/src/state.rs` | **Temperature is derived, never stored.** You cannot set it; you add energy and the substance decides. Sometimes it decides nothing at all, and melts instead. |
| `crates/au-physics/src/grid.rs` | Finite-volume lattice. Cells hold *extensive* totals (you cannot conserve a density). |
| `crates/au-physics/src/transport.rs` | Conduction (flux-form, harmonic-mean face conductivity), Stefan–Boltzmann radiation, hydrostatic pressure, CFL stability + substepping. |
| `crates/au-physics/tests/physics.rs` | 12 tests. Conservation + validation against Fourier, Stefan, Clausius–Clapeyron. |
| `crates/au-sim/src/physics_columns.rs` | Column ids `0x1xxx`. Permanent — a save names its columns by number. |
| `crates/au-sim/src/systems/physics.rs` | The `Layer::Physics` system. Owns the **active set** guard and the **clock cap**. |
| `crates/au-sim/tests/thermodynamics.rs` | 9 integration tests: conservation end-to-end, save/resume with thermal history, and the active-set guard. |
| `data/reference_materials.kv` | **Validation fixture, not content.** The engine never loads it. Real measured Earth substances, used to check the engine against three centuries of known answers. Phase 3 deletes it. |

### Modified files

| File | Change | Why |
|---|---|---|
| `Cargo.toml` | Added `au-physics` to the workspace | — |
| `data/constants.kv` | Physics section: grid, gravity, boundaries, CFL budget | Constants are data. "What if gravity were different?" must not need a recompile. |
| `crates/au-core/src/event.rs` | Added `PHYSICS_TIME_CAPPED` | The engine needs a way to *say* it could not do what was asked. |
| `crates/au-data/src/codec.rs` | `i128` support; `FORMAT_VERSION` 1 → 2 | Exact quantities are 16 bytes wide. |
| `crates/au-data/src/column.rs` | `Field` impls for `i128`, `u16` | Physics columns. |
| `crates/au-data/src/snapshot.rs` | Persist the **active set** and the **ledger** | Neither is derivable from the chunks. Lose the active set and a reloaded world quietly stops simulating half of itself. |
| `crates/au-data/tests/memory_strategy.rs` | Updated for the v2 `encode` signature | Mechanical. |
| `crates/au-sim/Cargo.toml` | Depends on `au-physics` | — |
| `crates/au-sim/src/lib.rs` | Exports `physics_columns` | — |
| `crates/au-sim/src/config.rs` | `merge()`, `as_map()` | Lets the validation fixture be layered on without the engine shipping a substance of its own. |
| `crates/au-sim/src/world.rs` | Added `materials`, `active`, `ledger`; `activate()` | See "the active set", below. |
| `crates/au-sim/src/sim.rs` | Registers the physics system; snapshot v2 | — |
| `crates/au-sim/src/systems/mod.rs` | Added `physics` | — |
| `crates/au-headless/Cargo.toml` | Depends on `au-physics` | — |
| `crates/au-headless/src/main.rs` | New `au physics` demo | A 640 m column of rock, heated from below, radiating to space. |
| `README.md` | Phase 2 | — |

### The active set — a bug caught before it was written

Phase 1 promised that *looking* at the world does not change it. Physics destroys that promise silently unless you are careful: if the physics system ran on "whatever chunks are in memory", then the set of things being **simulated** would be the set of things somebody **looked at**. The universe would evolve differently depending on where the camera was pointed. Observation would become interaction, and every saved world would diverge from its own replay.

So `World::active` is an explicit, snapshotted, hashed set. Physics runs on those chunks and only those. `World::chunk()` (looking) does not enrol; `World::activate()` (simulating) does. Test: `observing_the_world_does_not_simulate_it` sweeps a camera across fifty chunks of real rock in a running world and the history is bit-identical.

A consequence worth stating: an active chunk is **dirty forever**. Physics wrote to it, so it is no longer a function of the seed and can never be evicted. That is not a flaw in the Phase 1 memory strategy — it is the memory strategy telling the truth about what an evolving world costs. It is also exactly why LOD is not optional.

### The finding that Phase 2 forced on Phase 1

**You cannot fast-forward a diffusion equation.**

An explicit heat solver is stable only while `dt ≤ dx²/(2·ndim·D)`. Past that it does not lose accuracy — it oscillates to infinity in about thirty steps. Doubling the seconds-per-tick does not make heat spread twice as fast.

So the physics layer now **caps the clock** at what it can actually integrate, and logs `PHYSICS_TIME_CAPPED`. The deep-time accelerator keeps pushing; physics keeps refusing; the whole thing plateaus at the truth. The engine reports that it cannot resolve the physics rather than quietly producing a number.

The uncomfortable consequence: **the deep-time accelerator built in Phase 1 is largely inert on a world with real physics in it.** Reaching evolutionary timescales will need *coarse-grained* physics at low LOD — an implicit solver, or a statistical model — not the same physics run recklessly. That is now a known Phase 2b/4 task rather than a Phase 8 catastrophe. It also validates the concern I flagged when the accelerator was built.

### Bugs found by the tests

**The engine knew more physics than I did.** The latent-heat test failed: I asserted the melt plateau would sit at water's textbook 273.15 K. The engine produced 273.1575 K. It was right and I was wrong — the test box is a *vacuum*, `melt_k` is quoted at one atmosphere, and water's Clausius–Clapeyron slope is negative, so removing 101 kPa pushes its melting point 7.5 mK **up**. The engine had already applied that, correctly, to four decimal places, while the test insisted on the melting point of ice in a room.

### Validation — the engine reproduces known physics

You cannot assert that something *interesting* happens; that is the whole point of the project. But physics has right answers, worked out over three centuries, and if the engine does not reproduce them then every emergent thing built on top is emerging from a lie.

* **Fourier (1822)** — steady-state conduction gives a linear profile of slope q/k. ✓
* **Stefan (1879)** — a radiating surface settles at `T = (q/εσ)^¼`. ✓
* **Latent heat** — the plateau lasts *exactly* m·L/P seconds. ✓
* **Clausius–Clapeyron** — pressure raises rock's melting point and *lowers* water's. ✓
* **Second law** — heat never flows uphill; no cell may spontaneously invert. A Maxwell's demon in Phase 2 is a perpetual-motion organism in Phase 9. ✓

### What it looks like

`au physics` — a 640 m column of rock, 10 W/m² in at the floor, radiating to 2.7 K at the sky, starting uniformly solid at 300 K. Nothing tells it there should be a molten layer.

```
t =       2649 yr   floor  1318.6 K   sky   91.2 K   all solid
t =       8831 yr   floor  2121.2 K   sky  108.8 K   MOLTEN to 115 m
t =      26493 yr   floor  2847.2 K   sky  116.3 K   MOLTEN to 225 m
t =      52986 yr   floor  3023.7 K   sky  118.1 K   MOLTEN to 245 m

CONSERVATION ERROR                     0 µJ   ← not "small". zero.
```

The solid region has a gradient of exactly 3.33 K/m. The molten region has 6.65 K/m — because liquid rock conducts half as well as solid rock, so the same heat needs twice the gradient to get through. **Nobody wrote that bend.** It appeared where the rock melted, which is where the thermodynamics put it.

The surface settles at 118.1 K. Stefan–Boltzmann, computed by hand, says 118.3 K.

### Not done, on purpose

**No momentum, no fluid advection, no convection.** Convection is where planets get interesting (mantle overturn → tectonics → volcanism → hydrothermal vents), and I want it. But momentum needs the same exact-conservation treatment energy just got, and building it on an unproven energy substrate would be reckless. It is Phase 2b, and it is a real omission rather than an oversight.

**No cross-chunk halo exchange.** Physics operates within a chunk; chunk edges are boundary conditions. The two-pass gather/scatter design already accommodates it — it becomes necessary in Phase 4, when a planet spans chunks.

---

## Phase 1 — Simulation Foundation (initial commit)

Roadmap Phase 1 asks for: project structure, universal simulation clock,
save/load, basic data architecture, procedural world container. All five are
here, plus one addition argued for below.

### Added

| File | What it is | Why |
|---|---|---|
| `Cargo.toml` | Workspace, 4 crates | Zero third-party deps in the sim core: every crate is a place a version bump can silently change an iteration order or a float path |
| `data/constants.kv` | The dials of the universe | Constants are data, not code — "what if the constants were different?" is a question this engine must answer without a recompile |
| **au-core** | | |
| `time.rs` | Exact 128-bit fixed-point clock | `f64` accumulation drifts; drift is divergence. Also: **simulated** seconds-per-tick (sim state) is separated from the player's fast-forward (presentation), so watching faster cannot alter history |
| `rng.rs` | Splittable derived RNG streams | No global RNG, ever. `Rng = f(seed, domain, where, when)` — nothing shared, so parallelism is free and determinism survives it |
| `layer.rs` | The 14-layer stack as an ordered type | Turns ARCHITECTURE's dependency chain from prose into something checkable |
| `schedule.rs` | Multi-rate scheduler + layer enforcement | ARCHITECTURE §3 (systems tick at different rates), §20 (temporal LOD), §22 (Golden Rule, as a boot-time assertion) |
| `ids.rs` | Index+generation entity handles | Family trees are a promised feature; a stale handle must be detectably stale, not silently point at whoever reused the slot |
| `event.rs` | Append-only causal record | VISION's Fundamental Rule is only worth something if we can *check* it |
| `hash.rs` | 64-bit world digest | The single most valuable diagnostic in the project. Written before there is anything to hash, because by Phase 8 the causal chain is 10⁹ ticks long and retrofitting is impossible |
| **au-data** | | |
| `column.rs` | SoA storage, type-erased | AoS→SoA is 1–2 orders of magnitude at these entity counts. Decides whether Phase 9 is a simulation or a slideshow |
| `chunk.rs` | Chunk store, LOD, dirty tracking | ARCHITECTURE §19/§20. Untouched world isn't data — it's a function of the seed. Only *changes* cost memory |
| `codec.rs` | Hand-rolled versioned binary format | Every data structure here will change shape before Phase 13. If old saves can't migrate, we start avoiding schema changes to protect our saves — that's how a codebase ossifies |
| `snapshot.rs` | Seed + dirty chunks + event log | 2,000,000 ticks save in 160 bytes |
| **au-sim** | | |
| `world.rs` | The substrate | Should stay roughly this size forever — layers arrive as *columns* and *systems*, not as fields here |
| `sim.rs` | The tick loop | Small on purpose. If it ever grows a `match` on layer, the architecture has failed and this is where it shows |
| `config.rs` | Tiny `key = value` parser | Same reasoning as zero-deps: this file's contents change the outcome of the universe |
| `systems/time_scale.rs` | Deep-time acceleration | Chemistry wants seconds, evolution wants millennia; the gap is ~10¹⁴ and no machine runs 10¹⁴ ticks |
| **au-headless** | | |
| `main.rs` | `run` / `verify` / `resume` / `bench` | You cannot assert "a wolf appears". The only assertions available are structural, and this is what checks them |

### Added beyond the roadmap

**Observability.** The roadmap doesn't list it; the project cannot survive
without it. You cannot debug emergence by looking at it — when Phase 8 produces
no life after 400 million simulated years, the causal chain is 10⁹ ticks long
and the *only* tool that helps is deterministic replay. So the world hash, the
event log and the `verify`/`resume` commands are Phase 1 deliverables, not
Phase 14 conveniences.

### Bugs found and fixed by the tests

1. **Eviction was emitting an event.** Events feed the world hash, so memory
   pressure would have become part of the identity of the universe — a host with
   less RAM would have computed a different history. Would have surfaced years
   later as "evolution runs slightly differently on the server."
   *Invariant now enforced: nothing the engine does to manage memory may be
   observable to the simulation.*

2. **The snapshot wrote the event log's `len` where its `cap` belongs.** A world
   saved with one event in the log came back with a log that could hold exactly
   one event — and began silently discarding history on the very next emit. A
   save file that quietly truncates the past is the one failure this project
   cannot survive. Caught by `a_saved_universe_resumes_into_the_same_future`.

### Tests: 31 passing

- **au-core (12)** — time doesn't drift over 10⁶ ticks; deep time reaches 4.5 Gyr;
  RNG streams are pure and mutually independent; a system that reads upward is
  *rejected at registration*; execution order is bottom-up regardless of
  registration order; recycled entity slots don't resurrect the dead.
- **au-data (9)** — evicted chunks regenerate bit-identically; modified chunks are
  never evicted; the world hash ignores memory pressure; floats survive the codec
  bit-exactly; an unregistered column is an error, not a shrug; a corrupted save
  refuses to load.
- **au-sim (10)** — same seed → same universe, *tick by tick*; a saved universe
  resumes into the same future; two hosts with opposite eviction policies compute
  the same history; the time scale survives a save; and
  `the_world_contains_nothing_because_nothing_exists_yet`, which fails the day
  someone gets impatient and seeds the world with placeholder scenery.

### Not done, on purpose

No physics, matter, chemistry, terrain, or life. No stub modules for them either.
Matter arrives in Phase 3; anything in the world before then is content, and
content is the one thing the Golden Rule forbids.

---

## Phase 2 — Thermodynamics

Roadmap Phase 2 asks for: time progression, gravity, temperature, energy,
materials. All five, plus the conservation machinery that makes them trustworthy.

### The decision that shaped everything else

**Physics here is thermodynamics, not mechanics.** The instinct is to build rigid
bodies and collisions. That would be the wrong thing built beautifully. Chemistry
needs temperature and available energy; climate *is* heat transport; and life is a
dissipative structure that exists because a free-energy gradient exists. Nobody in
this stack needs a falling crate.

### Added

| File | What | Why |
|---|---|---|
| `crates/au-physics/quantity.rs` | Exact `Energy(i128 µJ)`, `Mass(i128 µg)`, `Ledger` | Evolution is an adversarial optimiser. A float leak is a free-energy exploit, and something will evolve to eat it |
| `crates/au-physics/material.rs` | Bulk material model, phase transitions, Clausius–Clapeyron | Phase 2 says what a substance *does*; Phase 3 will say what it *is*. Ships zero substances |
| `crates/au-physics/state.rs` | Temperature as a **derived** quantity, with latent heat | You cannot set the temperature. You add energy. Storing T would run cause and effect backwards |
| `crates/au-physics/grid.rs` | Finite-volume lattice, boundary conditions | Cells hold *extensive* totals, not densities — you cannot conserve a density |
| `crates/au-physics/transport.rs` | Conduction, radiation, hydrostatics, CFL substepping | Flux form: one integer per face, one sign each way. Conservation by construction |
| `crates/au-sim/physics_columns.rs` | Column ids `0x1xxx` | Ids are permanent — a renumbered column loads *wrong*, not not-at-all |
| `crates/au-sim/systems/physics.rs` | The physics system + the active-set guard | See below |
| `data/reference_materials.kv` | Real Earth substances | **A measuring instrument, not content.** The engine never loads it |

### Modified

`Cargo.toml` (workspace) · `au-core/event.rs` (`PHYSICS_TIME_CAPPED`) ·
`au-data/codec.rs` (i128/u16, FORMAT_VERSION 1→2) · `au-data/column.rs` (`Field` for
i128/u16) · `au-data/snapshot.rs` (v2: active set + ledger) · `au-sim/world.rs`
(materials, active set, ledger) · `au-sim/sim.rs` · `au-sim/config.rs` (`merge`) ·
`au-headless/main.rs` (`physics` demo) · `data/constants.kv`

All ten Phase 1 tests pass **unchanged**. An empty world imposes no speed limit, so
nothing about Phase 1's behaviour moved.

### The bug the active set exists to prevent

Physics could have silently destroyed Phase 1's promise that *looking at the world
does not change it.* If the physics system ran on "whatever chunks are resident,"
then the set of things being **simulated** would be the set of things somebody
**looked at** — and the universe would evolve differently depending on where the
camera pointed. Observation would become interaction, and every saved world would
diverge from its own replay.

So `World::active` is explicit, snapshotted and hashed. `World::activate()` is the
only door in, and generating a chunk to look at it does not knock.

### The discovery

**You cannot fast-forward a diffusion equation.** Explicit heat solvers diverge past
`dt = dx²/(2·ndim·D)` — not degrade, diverge. So Phase 1's deep-time accelerator is
largely inert once there is real physics in the world. The physics layer now caps
the clock at what it can integrate and logs it.

Deep time will need *coarse-grained* physics at low LOD, not the same physics run
faster. That is a real constraint on the project. Better found in Phase 2 than
Phase 8.

Related: **an active chunk can never be evicted again** — physics wrote to it, so it
is no longer a function of the seed. An evolving world costs memory forever, which
is precisely why LOD is mandatory.

### The test that failed because the engine was right

`latent_heat_produces_a_plateau_of_exactly_the_right_length` compared against water's
textbook 273.15 K and found nothing. The engine had produced 273.1575 K — because
the test box is a *vacuum*, and water's Clausius–Clapeyron slope is negative, so
removing one atmosphere pushes its melting point 7.5 mK *up*. The engine had already
applied that, correctly, to four decimals. The test hadn't.

### Tests: 52 passing (was 43)

- **au-physics (12)** — energy conserved *exactly* over 128k substeps with heat
  flowing in and out; a sealed box changes by zero; Fourier's linear profile;
  Stefan–Boltzmann equilibrium (engine 118.1 K vs analytic 118.3 K); a latent-heat
  plateau of exactly m·L/P seconds; pressure moving the melting point, with the
  right *sign* for both rock and water; the second law; the engine admitting when it
  cannot resolve the timestep.
- **au-sim/thermodynamics (9)** — the books balance to the microjoule through
  chunks, columns, snapshots and the scheduler; **observing the world does not
  simulate it**; an evolving chunk survives eviction; a thermal history survives a
  save; physics refuses to let the clock outrun it; and
  `a_gradient_produces_structure_that_nobody_designed`.

### Not done, on purpose

No momentum, no advection, no buoyancy — therefore no convection. It matters
(mantle → tectonics → vents → life), and it is deferred deliberately: mass and
momentum need the same exact-conservation guarantees energy now has, and building
them on an unproven energy substrate would have been reckless.
