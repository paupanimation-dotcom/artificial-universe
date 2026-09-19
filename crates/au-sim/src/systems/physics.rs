//! The physics system.
//!
//! # The bug this file exists to prevent
//!
//! Phase 1 established that *looking* at the world must not change it — chunks
//! are generated lazily, and a chunk generated because the camera swung past is
//! not supposed to alter history. There is a test called
//! `touching_chunks_does_not_change_history`.
//!
//! Physics breaks that, immediately and silently, unless we are careful. If the
//! physics system simply ran on "every chunk that happens to be in memory", then
//! the set of chunks being *simulated* would be the set of chunks somebody
//! *looked at* — and the universe would evolve differently depending on where
//! the player pointed the camera. Observation would become interaction. Every
//! saved world would diverge from its own replay.
//!
//! So there is an **active set**: an explicit, snapshotted, hashed list of the
//! chunks under simulation. Physics runs on those and only those.
//! [`World::activate`] is the only way in, and generating a chunk for viewing
//! does not call it. In Phase 4, the planet will decide what is active. Here,
//! the demo does.
//!
//! # The second thing this file does: it stops the clock lying
//!
//! Diffusion has a speed limit (see `au_physics::transport`). When the deep-time
//! accelerator asks for a timestep that cannot be integrated, this system
//! **slows the clock down to what the physics can actually do** and records the
//! fact. The accelerator keeps pushing, physics keeps refusing, and the whole
//! thing plateaus at the truth.
//!
//! That is a genuine finding, and it is inconvenient: it means the deep-time
//! accelerator built in Phase 1 is largely inert on a world with real physics in
//! it. Reaching evolutionary timescales will require *coarse-grained* physics —
//! an implicit solver, or a statistical model at low LOD — not the same physics
//! run faster. Better to know now.

use au_core::event::kind;
use au_core::layer::Layer;
use au_core::schedule::{System, SystemDesc};
use au_core::time::SimDuration;
use au_data::chunk::ChunkCoord;
use au_physics::{
    step, Boundary, Energy, FluidView, GridSpec, GridView, Momentum, Scratch, Wall,
};

use crate::physics_columns::{
    has_fluid, has_physics, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL, PHYS_MOM_X, PHYS_MOM_Y,
    PHYS_MOM_Z, PHYS_PRESSURE,
};
use crate::world::World;

pub struct PhysicsSystem {
    spec: GridSpec,
    boundary: Boundary,
    max_substeps: u32,
    safety: f64,

    // Working memory. All of it derived, none of it state. If clearing any of
    // this between ticks changed the result, determinism would already be gone.
    //
    // Gather → compute → scatter. Now that the fluid moves matter and momentum
    // around, seven columns have to be handed to the solver at once, and the
    // borrow checker will not let seven slices out of one `ColumnSet`. Copying
    // them out is not a workaround — it is the shape a SIMD or GPU kernel would
    // want anyway, and it is where that will slot in.
    scratch: Scratch,
    mass: Vec<i128>,
    material: Vec<u16>,
    energy: Vec<i128>,
    mom: [Vec<i128>; 3],
    press: Vec<f64>,

    // Telemetry.
    pub last_substeps: u32,
    pub last_stable_dt: f64,
    pub capped: bool,
}

impl PhysicsSystem {
    pub fn from_config(world: &World) -> Option<(SystemDesc, PhysicsSystem)> {
        if !world.config.bool_or("physics.enabled", true) {
            return None;
        }
        let c = &world.config;
        let spec = GridSpec::new(
            c.u64_or("physics.grid.nx", 1) as usize,
            c.u64_or("physics.grid.ny", 1) as usize,
            c.u64_or("physics.grid.nz", 1) as usize,
            c.f64_or("physics.grid.cell_m", 10.0),
        );
        let wall = |k: &str, d: &str| {
            Wall::parse(c.str_or(k, d)).unwrap_or_else(|| panic!("bad wall type for {}", k))
        };
        let boundary = Boundary {
            bottom_flux: c.f64_or("physics.bottom_flux_w_m2", 0.0),
            top_radiates: c.bool_or("physics.top_radiates", true),
            space_k: c.f64_or("physics.space_temperature_k", 2.7),
            surface_pressure: c.f64_or("physics.surface_pressure_pa", 0.0),
            gravity: c.f64_or("physics.gravity_m_s2", 0.0),
            bottom_temp: c.f64_opt("physics.bottom_temp_k"),
            top_temp: c.f64_opt("physics.top_temp_k"),
            fluid: c.bool_or("physics.fluid.enabled", false),
            x_wall: wall("physics.fluid.x_wall", "noslip"),
            y_wall: wall("physics.fluid.y_wall", "noslip"),
            z_wall: wall("physics.fluid.z_wall", "noslip"),
            projection_iters: c.u64_or("physics.fluid.projection_iters", 60) as u32,
            sor_omega: c.f64_or("physics.fluid.sor_omega", 1.85),
            // Implicit conduction: off by default, so every earlier phase
            // reproduces bit-for-bit. Turning it on frees the timestep from the
            // dx² diffusion limit — the change that makes deep time reachable.
            implicit_conduction: c.bool_or("physics.implicit_conduction", false),
            conduction_iters: c.u64_or("physics.conduction_iters", 24) as u32,
            conduction_omega: c.f64_or("physics.conduction_omega", 1.0),
            viscous_iters: c.u64_or("physics.viscous_iters", 0) as u32,
        };

        let desc = SystemDesc::new("physics.transport", Layer::Physics)
            .reads(&[Layer::Universe, Layer::Physics])
            .every(1);

        Some((
            desc,
            PhysicsSystem {
                spec,
                boundary,
                max_substeps: c.u64_or("physics.max_substeps", 256) as u32,
                safety: c.f64_or("physics.cfl_safety", 0.4),
                scratch: Scratch::new(),
                mass: Vec::new(),
                material: Vec::new(),
                energy: Vec::new(),
                mom: [Vec::new(), Vec::new(), Vec::new()],
                press: Vec::new(),
                last_substeps: 0,
                last_stable_dt: f64::INFINITY,
                capped: false,
            },
        ))
    }

    pub fn spec(&self) -> GridSpec {
        self.spec
    }
    pub fn boundary(&self) -> Boundary {
        self.boundary
    }
}

impl System<World> for PhysicsSystem {
    fn run(&mut self, world: &mut World) {
        // An empty active set means there is nothing under simulation — which is
        // the shipped state of this engine, and the reason all of Phase 1's
        // tests still pass untouched.
        if world.active.is_empty() {
            return;
        }
        let dt = world.clock.scale().as_secs_f64();
        if dt <= 0.0 {
            return;
        }

        let coords: Vec<ChunkCoord> = world.active.iter().copied().collect();
        // Make sure every active chunk is resident. Note this is the *simulation*
        // asking for them, not a camera.
        for c in &coords {
            world.chunk(*c);
        }

        // Disjoint field borrows: chunks mutably, materials immutably. Events and
        // the ledger are accumulated locally and folded in afterwards, because
        // `World::emit` needs the whole world.
        let mats = &world.materials;
        let chunks = &mut world.chunks;

        let b_fluid = self.boundary.fluid;
        let mut impulse = [0i128; 3];
        let mut e_in: i128 = 0;
        let mut e_out: i128 = 0;
        let mut shortest_dt = f64::INFINITY;
        let mut unresolved = false;
        let mut substeps = 0u32;
        let mut stable_dt = f64::INFINITY;

        for c in &coords {
            let Some(chunk) = chunks.get_mut(*c) else { continue };
            if !has_physics(chunk.columns()) {
                continue;
            }

            // ── Gather. ─────────────────────────────────────────────────────
            let fluid_here = b_fluid && has_fluid(chunk.columns());
            {
                let cols = chunk.columns();
                self.mass.clear();
                self.mass.extend_from_slice(cols.get::<i128>(PHYS_MASS).unwrap().as_slice());
                self.material.clear();
                self.material
                    .extend_from_slice(cols.get::<u16>(PHYS_MATERIAL).unwrap().as_slice());
                self.energy.clear();
                self.energy
                    .extend_from_slice(cols.get::<i128>(PHYS_ENERGY).unwrap().as_slice());
                if fluid_here {
                    for (k, id) in [PHYS_MOM_X, PHYS_MOM_Y, PHYS_MOM_Z].iter().enumerate() {
                        self.mom[k].clear();
                        self.mom[k]
                            .extend_from_slice(cols.get::<i128>(*id).unwrap().as_slice());
                    }
                    self.press.clear();
                    self.press
                        .extend_from_slice(cols.get::<f64>(PHYS_PRESSURE).unwrap().as_slice());
                }
            }
            if self.mass.len() != self.spec.cells() {
                // Built for a different grid. Refuse rather than integrate garbage.
                continue;
            }

            // ── Compute. ────────────────────────────────────────────────────
            let mut boundary = self.boundary;
            boundary.fluid = fluid_here;

            let rep = {
                let (m0, m1, m2) = {
                    let [a, b, c] = &mut self.mom;
                    (a, b, c)
                };
                let mut view = GridView::new(
                    self.spec,
                    &mut self.mass,
                    &self.material,
                    &mut self.energy,
                );
                if fluid_here {
                    view = view.with_fluid(FluidView {
                        mom_x: m0,
                        mom_y: m1,
                        mom_z: m2,
                        pressure: &mut self.press,
                    });
                }
                step(
                    &mut view,
                    mats,
                    &boundary,
                    dt,
                    self.max_substeps,
                    self.safety,
                    &mut self.scratch,
                )
            };

            // ── Scatter. `columns_mut` dirties the chunk — correctly. A chunk
            // physics has touched is no longer a function of the seed and can
            // never be evicted again. Not a flaw in the memory strategy; the
            // memory strategy telling the truth about what an evolving world
            // costs. It is also exactly why LOD is not optional.
            {
                let cols = chunk.columns_mut();
                cols.get_mut::<i128>(PHYS_MASS)
                    .unwrap()
                    .as_mut_slice()
                    .copy_from_slice(&self.mass);
                cols.get_mut::<i128>(PHYS_ENERGY)
                    .unwrap()
                    .as_mut_slice()
                    .copy_from_slice(&self.energy);
                if fluid_here {
                    for (k, id) in [PHYS_MOM_X, PHYS_MOM_Y, PHYS_MOM_Z].iter().enumerate() {
                        cols.get_mut::<i128>(*id)
                            .unwrap()
                            .as_mut_slice()
                            .copy_from_slice(&self.mom[k]);
                    }
                    cols.get_mut::<f64>(PHYS_PRESSURE)
                        .unwrap()
                        .as_mut_slice()
                        .copy_from_slice(&self.press);
                }
            }

            for a in 0..3 {
                impulse[a] += rep.impulse[a];
            }
            e_in += rep.energy_in.0;
            e_out += rep.energy_out.0;
            substeps = substeps.max(rep.substeps);
            stable_dt = stable_dt.min(rep.stable_dt);
            if rep.unresolved {
                unresolved = true;
                shortest_dt = shortest_dt.min(rep.dt_done);
            }
        }

        world.ledger.energy_in += Energy(e_in);
        world.ledger.energy_out += Energy(e_out);
        for a in 0..3 {
            world.ledger.impulse[a] += Momentum(impulse[a]);
        }
        self.last_substeps = substeps;
        self.last_stable_dt = stable_dt;
        self.capped = unresolved;

        if unresolved && shortest_dt.is_finite() && shortest_dt > 0.0 {
            // The physics could not advance as far as the clock wanted. Rather
            // than let the clock run ahead of the world it describes — which
            // would make every timestamp in the history a lie — we slow the
            // clock to what actually happened.
            //
            // The deep-time accelerator will try again next period, and be
            // refused again. It plateaus here, at the truth, and the event log
            // says so.
            let capped = SimDuration::from_secs_f64(shortest_dt);
            world.clock.set_scale(capped);
            world.emit(
                Layer::Physics,
                kind::PHYSICS_TIME_CAPPED,
                capped.secs,
                capped.frac,
                substeps as u64,
            );
        }
    }
}
