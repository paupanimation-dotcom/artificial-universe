# Artificial Universe

> The engine should never ask *"What content should I create?"*
> It should ask *"What happens if these rules interact for enough time?"*
> — ARCHITECTURE §22

**Status: Phase 5i complete — protocells divide because of what they are made of. 237 tests passing. Conservation error: exactly 0.**

---

## Quick start

```bash
cargo test --release          # 237 tests
cargo run --release -- physics --ticks 3000   # a slab of rock, a heat source, a cold sky
cargo run --release -- convect                # Rayleigh–Bénard rolls from a heated fluid
cargo run --release -- react                  # the au-chem reactor directly: ignition + equilibrium
cargo run --release -- react-sim              # chemistry IN the sim loop: shared energy field, save/resume
cargo run --release -- planet                 # a world from a star, a mass, an orbit; the habitable zone, emergent
cargo run --release -- deeptime               # implicit conduction: 5000× fewer substeps, still 0 µJ error
cargo run --release -- genesis                # H and O atoms in a hot box; the engine discovers water
cargo run --release -- vent                   # water invented at one cell, found at the tube's ends
cargo run --release -- raf                    # catalysis, and a molecule that makes more of itself
cargo run --release -- verify  --ticks 500000 # same seed → same universe?
cargo run --release -- resume  --ticks 500000 --at 180000
```

## The stack so far

```
crates/
  au-core/      time · RNG streams · layers · scheduler · events · hashing
  au-data/      SoA columns · chunk store · LOD · codec · snapshots
  au-physics/   thermodynamics · fluid dynamics · implicit conduction · species diffusion · exact conservation
  au-chem/      elements (conserved) · molecules as graphs · generative reaction network · reactor
                · physical rate constants · derived catalysis · RAF detection · amphiphilicity
                · self-assembly · open boundaries and the atom ledger · vesicle geometry
  au-planet/    planets as parameters · derived state · seeded relief (documented)
  au-sim/       World · config · tick loop · the active set
  au-headless/  CLI: run, verify, resume, bench, physics, convect, ra, react, react-sim, planet, deeptime, genesis, vent, raf
data/
  constants.kv              the dials of the universe — data, not code
  reference_materials.kv    a MEASURING INSTRUMENT, not content (see below)
```

**Zero third-party dependencies in the simulation core.** Every crate is a place
where a version bump can silently change an iteration order or a floating-point
path, and determinism is what this whole project rests on.

---

## Phase 1 — the deterministic foundation

**Time is exact and is not a float.** `u64` seconds + `u64` fraction (2⁻⁶⁴ s).
Adding `0.1` a million times in `f64` gives you `100000.00000133288`. That drift is
a divergence, and divergence is death. Fixed point cannot drift.

**There is no global RNG.** Randomness is *derived*: `Rng = f(seed, domain, where,
when)`. Nothing shared ⇒ ordering is irrelevant ⇒ parallelism is free and
determinism survives it.

**The layer stack is a boot-time assertion.** Every system declares what it writes
and what it reads. Read *upward* — physics consulting ecology — and the engine
refuses to start.

**Untouched world is not data — it is a function of the seed.** Only what the
simulation *changed* costs memory or disk.

## Phase 2 — thermodynamics

### What "physics" means here

Not rigid bodies. Not collisions. Look at what the layers above actually need:
chemistry needs temperature and available energy; climate *is* heat transport; and
the origin of life needs one thing above all — a **persistent free-energy
gradient**, because life is a dissipative structure. It exists because the star is
hot and the sky is cold. It is what grows in the gap.

Nobody in this stack needs a falling crate. Everybody needs *how much energy is
here, and where is it going.* Building a rigid-body engine would have been building
the wrong thing beautifully.

### Conservation is exact, and that is not fussiness

**Evolution is an adversarial optimiser whose search space includes the bugs in
your physics.** Leave any path by which energy can be created and something will
eventually evolve to walk it — not because it is clever, but because natural
selection is an exhaustive search and free energy is the most rewarding thing in
any possible fitness landscape. This is how artificial-life projects die.

So energy and mass are **exact integers** (µJ and µg, `i128`) and transport is
**symmetric flux** — the same integer subtracted from one cell and added to the
other. Conservation is an identity, not an aspiration:

```
E(now) − E(start)  ==  energy_in − energy_out       exactly, forever
```

Floats still appear — temperature, density, pressure — but only as **derived**
values, recomputed from the exact integers every tick and never stored. *Float error
cannot accumulate if you never accumulate floats.*

### You cannot set the temperature

You add energy, and the temperature is what you *get*. Sometimes — and this is the
good part — the substance decides on nothing at all, and melts instead.

That plateau is a thermostat with no thermostat in it. It is why Earth has a
climate rather than a temperature swing, and it is the first emergent regulation in
this codebase: nobody wrote a stabiliser.

---

## The demo

```
$ cargo run --release -- physics --ticks 3000

  a 640 m column of rock, 64 cells
  10 W/m² in at the floor, radiating to 2.7 K at the sky
  starts uniformly solid at 300 K

  Nothing tells it there should be a molten layer.

  t =          0 yr   floor   300.0 K   sky  300.0 K   all solid
  t =       2649 yr   floor  1318.6 K   sky   91.2 K   all solid
  t =       8831 yr   floor  2121.2 K   sky  108.8 K   MOLTEN to 115 m
  t =      26493 yr   floor  2847.2 K   sky  116.3 K   MOLTEN to 225 m
  t =      52986 yr   floor  3023.7 K   sky  118.1 K   MOLTEN to 245 m

  CONSERVATION ERROR                     0 µJ   ← not 'small'. zero.
```

Look at the kink. The solid region has a gradient of exactly **3.33 K/m**; the
molten region has **6.65 K/m** — because liquid rock conducts half as well as solid
rock, so the same heat needs twice the gradient to get through.

Nobody wrote that bend. It appeared where the rock melted, which is where the
thermodynamics put it. That is the entire thesis of the project, at the smallest
scale it can be demonstrated: **put in laws, get out structure.**

---

## The three invariants (Phase 1) and the fourth (Phase 2)

| Invariant | Why it decides everything |
|---|---|
| Same seed → same universe, tick by tick | When Phase 8 produces no life after 400M years, replay is the only way to find out why |
| A saved universe resumes into the same future | Otherwise every family tree the player inspects is a fabrication |
| Memory management is invisible to the world | Otherwise a host with less RAM computes a different history |
| **Observation is not interaction** | Otherwise the universe evolves differently depending on where the camera points |

That fourth one is new and it is why `World::active` exists. Physics runs on the
**active set** — an explicit, snapshotted, hashed list of chunks under simulation —
not on "whatever happens to be in memory." Generating a chunk to *look* at it does
not enrol it.

## What Phase 2 validated against reality

The engine is handed real, measured substances and asked to reproduce results known
before computers existed. `data/reference_materials.kv` is a **measuring
instrument**, not content — the engine never loads it, and Phase 3 will delete it
once materials are *derived* from elements.

- linear conduction profile — **Fourier, 1822** ✓
- surface at T = (q/εσ)^¼ — **Stefan, 1879**: engine 118.1 K, analytic 118.3 K ✓
- latent-heat plateau lasting exactly m·L/P seconds ✓
- melting point moving with pressure, and moving the *wrong way* for water ✓
- the second law: heat never flows uphill, no Maxwell's demons ✓

That last water test failed at first — **because the engine knew more physics than
I did.** I compared against 273.15 K, the textbook number. The engine produced
273.1575 K, because the box is a vacuum, and water's Clausius–Clapeyron slope is
negative, so removing an atmosphere pushes its melting point *up* by 7.5 mK. It had
already applied that. My test hadn't.

---

## The uncomfortable discovery

**You cannot fast-forward a diffusion equation.**

An explicit heat solver is stable only while `dt ≤ dx²/(2·ndim·D)`. Past that it
does not lose accuracy — it oscillates to infinity in about thirty steps. Doubling
the seconds-per-tick does not make heat spread twice as fast.

Which means the deep-time accelerator built in Phase 1 is **largely inert on a
world with real physics in it.** The physics layer now caps the clock at what it can
actually integrate and logs the fact (`PHYSICS_TIME_CAPPED`). The accelerator keeps
pushing; physics keeps refusing; the whole thing plateaus at the truth.

Reaching evolutionary timescales will require *coarse-grained* physics at low LOD —
an implicit solver, or a statistical model — not the same physics run recklessly.
That is a real constraint on the whole project, and it is now a known one rather
than a surprise waiting in Phase 8.

A second, related cost surfaced with it: **an active chunk can never be evicted
again.** Physics wrote to it, so it is no longer a function of the seed. An evolving
world costs memory forever — which is exactly why LOD is not optional. Most of a
planet must *not* be finely evolving, or the planet will not fit in the machine.

---

## Open questions for Phase 3+

- **Chemistry cannot be literal.** Real abiogenesis is not computable at any scale
  we can run. Phases 3–5 need an abstraction that preserves causal structure while
  staying tractable — a reaction-network model, not molecular dynamics. Choosing it
  consciously is the main design work ahead.
- **No momentum yet.** Nothing moves matter: no advection, no buoyancy, therefore no
  convection. Convection matters enormously (mantle → tectonics → volcanoes →
  hydrothermal vents), and when it arrives, mass and momentum get the same
  exact-integer symmetric-flux treatment energy already has. Deferred deliberately:
  building it on an unproven energy substrate would have been reckless.
- **LOD transitions must conserve quantities.** If demoting a chunk loses biomass,
  ecosystems will drift every time the camera moves. Every transition must be a
  *summary*, never a truncation.

## Phase 5i step 2: composition determines replication rate

A compartment that merely persists is a rock. Phase 5h showed a membrane is
worth having in a flow; step 1 showed when a bag has enough surface to become
two. Neither gives anything to inherit. **Without a link from what a bag is made
of to how fast it divides, composition is inherited information with no
consequences** — variation cannot affect persistence, nothing can adapt, and you
have a crystal with a membrane. That is the failure artificial life hits most
often, and it is the one thing this step had to beat.

```
events=335   surfactant inside a bag: 12,567 -> 12,902
v=0.6605  area=5.295e-15     <- below 1/sqrt(2), fission
v=0.9030  area=2.676e-15     <- daughter: area halved, a stable sphere again
```

A bag imports substrate, converts it to surfactant *inside itself*, and the
surfactant cannot leave — because it **is** the membrane. Area ratchets upward,
`v` crosses the two-sphere bound, the bag divides, and both daughters inherit
the composition that caused it. A world whose chemistry cannot build membrane
never divides at all.

Nothing in the code grants membranes an advantage. The advantage is what
trapping a product does inside a compartment. And the daughter landing at
v ≈ 0.90 is the √2 jump step 1 derived from geometry — observed, not asserted.

Protocells are individuals: id, parent, birth tick, sparse contents. Division
frees the parent's id and allocates two fresh ones, each recording the parent,
so lineage is a clean binary tree and no id ever means two things. It is the
first use of the `EntityId` generation counter Phase 1 added with the note that
a family tree must never point at a stranger.

**Species populations now live in two places** — grid columns hold the bulk,
protocells hold what they have taken up — and every accounting must sum both.

Selection is *enabled*, not observed. One bag dividing faster than another is
the precondition. Watching a composition spread through a population is 5j.

## Phase 5i: a bag divides when geometry says it can

A compartment that merely persists is a rock — no descendants, so nothing about
it can be inherited, so a good barrier and a bad one never accumulate any
difference between them.

The tempting rule is *divide when the contents exceed N*, and that number would
instantly become the most load-bearing constant in the engine. So the criterion
is taken from geometry instead. A bag's **area** comes from how many
amphiphiles it has assembled. Its **volume** comes from osmosis — a membrane
passes water and not solute, so it swells until it is as dilute as its
surroundings. Their ratio is the reduced volume, and two facts about it are
pure arithmetic:

```
    v > 1        impossible — a sphere holds the most volume per unit area
                 there is, so the bag bursts (after ~3% areal strain)

    v <= 1/√2    two spheres of the same total volume fit in this much
                 membrane — division is geometrically available
```

`1/√2` is not tuned. It is what `√π/√(2π)` equals.

**So a protocell sits between two failures.** Make osmolytes faster than
membrane and it bursts; make membrane faster than osmolytes and it divides.
Whether a bag divides, bursts or persists is decided by the ratio of two rates
in its own chemistry — not by its size, and not by anything anyone wrote down.

### Heredity without a genome

Division is the engine's first *stochastic* physics — deliberately, and with
the determinism story intact, because randomness here is derived from
`(seed, domain, where, when)` and a resumed world draws identical numbers.
Stochastic is not the same thing as nondeterministic.

No template, no copying, no sequence. Just a bag splitting in two — and a
daughter resembles its parent **fifty times more closely** than a composition
built from the same abundances rearranged. The proportions survive because the
coin is fair per molecule.

And fidelity is not a setting: a species present in `n` copies is transmitted
with relative error `1/(2√n)`. **Copy number *is* heredity's fidelity**, which
is the same reason a real cell carries one genome rather than one molecule of
each protein.

The limit is recorded as a test rather than left to be discovered:
compositional inheritance *decays*. Each split adds its own √n noise and
nothing corrects it, because a composition cannot be checked — only re-drawn.
Fixing that needs something that can be compared against, which is a template,
which is Phase 6.

## Phase 5h: the world is open, and a barrier becomes worth having

Every phase before this one ran in a box. Energy could cross a boundary from
Phase 2 onward; matter never could. That is why the Phase 5c autocatalytic
survey could only claim its set *would* sustain itself **if fed** — nothing had
ever fed anything, and a box has exactly one destiny.

A port is a Dirichlet condition on population, the twin of the fixed-temperature
walls the fluid solver has used since 2b. A **source** holds its cell at a fixed
composition; a **sink** draws its cell to zero. Nothing declares a rate — the
port states a concentration and the gradient decides the throughput. And the
sink is *blind*: it removes every species without preference, because a drain
that spared some molecules and not others would be a fitness function hidden in
a boundary condition.

Then this happens, and nobody wrote it:

```
same tube · same tracer · same 20,000 ticks · one open port

    tracer surviving the drain    with a membrane    328,945
                                  without                720
                                                    ───────────
                                                       457×
```

Membranes lower permeability (Phase 5g). A drain removes what reaches it (Phase
5h). Neither fact is about advantage — and together they are one anyway. The
advantage is what a barrier *is*, in the presence of a flow.

It is not selection yet: no heredity, no individuals, nothing that can
differentially reproduce. It is the precondition — an environment where two
chemistries of identical composition have measurably different prospects,
decided by physics rather than by a parameter.

**Conservation does not weaken when the world opens.** `atoms(now) − atoms(start)
== in − out`, exactly, *per element* — bulk mass would balance perfectly through
a bug that turned carbon into oxygen. "Totals never change" was never the real
invariant; it was the special case where nothing crossed.

The validation anchor is Fick's steady solution in one dimension: a fixed source
and a fixed sink give a **linear** profile with the same flux across every face.
A 21-cell tube lands on the analytic ramp within 2%, holds its contents within
1% while continuing to pass matter through, and returns what it is fed within
2%. Which is the distinction the whole phase exists to draw: **equilibrium is
where a box ends; a steady state is somewhere a system can live.**

