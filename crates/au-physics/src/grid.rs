//! The grid.
//!
//! A finite-volume lattice. Cells hold *extensive* quantities — the total mass
//! and the total energy that happen to be inside that box — not intensive ones
//! like density or temperature. That choice is what makes conservation
//! expressible at all: you cannot conserve a density, but you can conserve a
//! total.
//!
//! `y` is up. Gravity points along −y. That is the only geometric convention in
//! the file and everything else follows from it.
//!
//! # No structs-of-cells
//!
//! There is no `Cell` type holding `{mass, energy, material}`. The data arrives
//! as three parallel slices, because that is how it lives in the chunk columns
//! (see `au-data`), and copying it into an array of structs every tick to make it
//! read nicer would cost more than the physics does.

use crate::material::MaterialId;

/// Dimensions and scale of a lattice. Cheap, `Copy`, no data.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct GridSpec {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    /// Edge length of a cell, metres. Cells are cubes.
    pub cell_m: f64,
}

impl GridSpec {
    pub fn new(nx: usize, ny: usize, nz: usize, cell_m: f64) -> Self {
        assert!(nx > 0 && ny > 0 && nz > 0, "grid must be non-degenerate");
        assert!(cell_m > 0.0 && cell_m.is_finite(), "cell size must be positive");
        GridSpec { nx, ny, nz, cell_m }
    }

    #[inline]
    pub fn cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Linear index. x fastest, then y, then z — so a vertical column is
    /// strided, and a horizontal row is contiguous. Heat mostly flows
    /// vertically in a planet, so this is arguably backwards; it is also
    /// arguably irrelevant until profiling says otherwise, and guessing at cache
    /// layouts before there is a workload is how you optimise the wrong thing.
    #[inline]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        debug_assert!(x < self.nx && y < self.ny && z < self.nz);
        x + self.nx * (y + self.ny * z)
    }

    /// Volume of one cell, m³.
    #[inline]
    pub fn cell_volume(&self) -> f64 {
        self.cell_m * self.cell_m * self.cell_m
    }

    /// Area of the face between two adjacent cells, m².
    #[inline]
    pub fn face_area(&self) -> f64 {
        self.cell_m * self.cell_m
    }

    /// How many axes actually have neighbours. A 1×64×1 column diffuses in one
    /// dimension, and the stability limit is three times more generous than it
    /// would be for a cube — so this is not a cosmetic detail.
    #[inline]
    pub fn dimensionality(&self) -> usize {
        (self.nx > 1) as usize + (self.ny > 1) as usize + (self.nz > 1) as usize
    }
}

/// What happens at the edge of the grid, per axis.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Wall {
    /// The grid wraps. A domain with no sides, which is how you study a fluid
    /// layer without the sides being the story.
    Periodic,
    /// Solid, and the fluid sticks to it. Tangential velocity → 0 at the surface.
    NoSlip,
    /// Solid, but frictionless. The fluid slides along it freely.
    ///
    /// Physically the odder of the two, and yet the one worth reaching for
    /// first: Rayleigh–Bénard convection between two stress-free surfaces has a
    /// **closed-form** critical Rayleigh number, 27π⁴/4 = 657.511. A solver that
    /// reproduces a number derived on paper in 1916 is a solver you can trust.
    /// No-slip has no closed form — only the tabulated 1707.762.
    FreeSlip,
}

impl Wall {
    pub fn parse(s: &str) -> Option<Wall> {
        match s {
            "periodic" => Some(Wall::Periodic),
            "noslip" | "no_slip" => Some(Wall::NoSlip),
            "freeslip" | "free_slip" => Some(Wall::FreeSlip),
            _ => None,
        }
    }
    #[inline]
    pub fn is_periodic(self) -> bool {
        self == Wall::Periodic
    }
}

/// Boundary conditions.
///
/// **This is a setup, not a law**, and the distinction matters. A real planet's
/// bottom heat flux comes from radioactive decay and primordial heat; Phase 4
/// will compute it from the planet's actual inventory. Until then it is a knob,
/// and it is a knob in exactly the way "the seed" is a knob: it is an initial
/// condition, not a scripted outcome.
///
/// The top radiates to space. That is not a knob. It is why a planet can have a
/// temperature at all — a body that only absorbed would heat without limit, and
/// a body that only radiated would freeze. Everything interesting lives in
/// between, in the steady state where in equals out. **Life is what grows in that
/// gap.**
#[derive(Clone, Copy, Debug)]
pub struct Boundary {
    /// Heat entering through the floor, W·m⁻².
    pub bottom_flux: f64,
    /// Does the top surface radiate to space?
    pub top_radiates: bool,
    /// Temperature of the sink, K. The cosmic microwave background is 2.7 K.
    pub space_k: f64,
    /// Pressure above the top cell, Pa. Zero for vacuum.
    pub surface_pressure: f64,
    /// Downward acceleration, m·s⁻². Phase 4 derives this from planetary mass;
    /// here it is given.
    pub gravity: f64,

    /// Fixed temperature at the floor / ceiling, K.
    ///
    /// A Dirichlet condition, and the one Rayleigh–Bénard is defined by: hold the
    /// bottom hot and the top cold and *do not tell the fluid what to do about
    /// it*.
    pub bottom_temp: Option<f64>,
    pub top_temp: Option<f64>,

    /// Is momentum live? When false the grid is a rigid conductor and this is
    /// Phase 2 exactly as it was.
    pub fluid: bool,

    pub x_wall: Wall,
    pub y_wall: Wall,
    pub z_wall: Wall,

    /// Fixed-iteration count for the pressure solve.
    ///
    /// Fixed, not tolerance-based, because a loop that stops "when converged"
    /// makes the number of iterations depend on the state, which makes the answer
    /// depend on the history in a way nothing else in this engine does. The
    /// residual is *reported* instead of being chased.
    pub projection_iters: u32,
    /// SOR over-relaxation factor. ~1.85 for grids of this size; 1.0 is plain
    /// Gauss–Seidel and converges far more slowly.
    /// Over-relaxation for the pressure projection.
    ///
    /// **Leave this at 0 and it is derived**, which is what it should be: the
    /// optimal factor for a Laplace problem is `2/(1+sin(π/n))`, a property of
    /// the *grid* rather than a preference. It rises toward 2 as the grid
    /// refines — 1.59 at n=12, 1.77 at n=24, 1.86 at n=41, 1.91 at n=64.
    ///
    /// The declared 1.85 this defaulted to is the optimum for a 41-cell grid,
    /// which is the tank Phase 2b validated Ra_c on, and it was silently wrong
    /// for every other grid since. Measured on a 12-cell grid at twenty sweeps:
    ///
    /// ```text
    ///     omega 1.00  ->  83.4% divergence
    ///     omega 1.40  ->  15.4%
    ///     omega 1.59  ->   1.8%     <- derived
    ///     omega 1.85  -> 185.1%     <- the declared constant
    ///     omega 1.95  -> 168.3%
    /// ```
    ///
    /// Past the optimum, SOR does not converge slowly. It oscillates — which is
    /// why 1.85 was *worse than no over-relaxation at all* on a grid it did not
    /// fit. A hundredfold error, in a constant nobody derived.
    ///
    /// A nonzero value is honoured, for tests that need to pin it.
    pub sor_omega: f64,

    /// Solve conduction implicitly (backward Euler) instead of explicitly.
    ///
    /// When true the `dx²` diffusion limit no longer constrains the timestep —
    /// the single change that makes deep time reachable. Accuracy still degrades
    /// as the step grows; stability does not. Off by default so every earlier
    /// phase reproduces bit-for-bit.
    pub implicit_conduction: bool,
    /// Fixed sweep count for the implicit conduction solve. Fixed rather than
    /// tolerance-based for the same reason as `projection_iters`.
    pub conduction_iters: u32,
    /// Relaxation factor for the conduction solve. **Defaults to 1.0 (plain
    /// Gauss–Seidel), and raising it trades away a guarantee.**
    ///
    /// The conduction update is a weighted average of a cell's own previous
    /// temperature and its neighbours', with positive weights that sum to one — a
    /// convex combination. That single fact means the iterate can never leave the
    /// range of temperatures it started from, *however few sweeps are run*: the
    /// discrete maximum principle holds even on a badly under-converged solve.
    ///
    /// Over-relaxation (ω > 1) converges faster but is no longer a convex
    /// combination, so an under-converged field may overshoot — and since the
    /// energy transported is that temperature difference multiplied by a very
    /// large dt, the overshoot is amplified rather than damped. It is exactly the
    /// failure mode this phase's stability test caught: energy still conserved to
    /// the microjoule, but heat piled up in places physics forbids.
    ///
    /// Faster convergence is not worth an unbounded field, so the default is safe.
    pub conduction_omega: f64,

    /// Gauss–Seidel sweeps for the implicit viscous solve. **0 disables it**,
    /// which is the default and reproduces every earlier phase byte for byte.
    ///
    /// Explicit viscosity is stable only while `dt ≤ dx²/(2·ndim·ν)`, and that
    /// `dx²` is what puts a micron-scale grid out of reach — the same cliff
    /// Phase 4b removed for heat. Setting this takes the viscous term off the
    /// timestep, at the cost of accuracy rather than stability.
    pub viscous_iters: u32,
}

impl Default for Boundary {
    fn default() -> Self {
        Boundary {
            bottom_flux: 0.0,
            top_radiates: false,
            space_k: 2.7,
            surface_pressure: 0.0,
            gravity: 0.0,
            bottom_temp: None,
            top_temp: None,
            fluid: false,
            implicit_conduction: false,
            conduction_iters: 24,
            conduction_omega: 1.0,
            viscous_iters: 0,
            x_wall: Wall::NoSlip,
            y_wall: Wall::NoSlip,
            z_wall: Wall::NoSlip,
            projection_iters: 60,
            sor_omega: 0.0, // 0 = derive from the grid
        }
    }
}

/// Optimal SOR factor for a Laplace problem on a grid whose smallest active
/// dimension is `n` cells: `2/(1+sin(π/n))`. Not tuned — this is the standard
/// result, and it is what `sor_omega = 0` asks for.
pub fn derived_sor_omega(n: usize) -> f64 {
    if n < 2 {
        return 1.0;
    }
    2.0 / (1.0 + (std::f64::consts::PI / n as f64).sin())
}

impl Boundary {
    #[inline]
    pub fn wall(&self, axis: usize) -> Wall {
        [self.x_wall, self.y_wall, self.z_wall][axis]
    }
}

/// The velocity field, on a **staggered (MAC) grid**.
///
/// # Why the velocities do not live where the cells do
///
/// The obvious layout — velocity and pressure both at cell centres — is broken,
/// and broken in a way that looks like physics. On a co-located grid the discrete
/// Laplacian in the pressure equation only ever couples cells *two* apart, so the
/// odd and even cells decouple completely. The solver then happily converges to a
/// pressure field with a checkerboard superimposed on it: alternating high/low,
/// invisible to the equations, and it drives spurious flow. People have published
/// convection cells that were checkerboard artefacts.
///
/// Harlow and Welch's 1965 fix is to put the velocities **on the faces**. Then
/// the divergence at a cell uses its own two faces, and the pressure gradient at
/// a face uses its own two cells, and the resulting Laplacian is compact. No
/// decoupling, no checkerboard, no shortcut.
///
/// # The index convention
///
/// `mom_x[i]` is the x-momentum on the **−x face of cell i** — the face between
/// cell (x−1, y, z) and cell (x, y, z). Same length as the cell arrays. The far
/// face of the last cell is not stored: on a periodic axis it *is* face 0, and on
/// a walled axis the velocity there is zero anyway.
pub struct FluidView<'a> {
    pub mom_x: &'a mut [i128],
    pub mom_y: &'a mut [i128],
    pub mom_z: &'a mut [i128],

    /// Dynamic pressure, Pa. **Solver state, not a fundamental field.**
    ///
    /// In the exact solution, pressure is a Lagrange multiplier determined
    /// entirely by the current velocity — carrying it between steps would be
    /// meaningless. With a *fixed-iteration* solve it is not exact, so the
    /// previous step's answer is a much better starting guess than zero, and the
    /// solve becomes several times cheaper. That makes it real state: it is
    /// snapshotted, it is hashed, and a resumed world gets the same warm start a
    /// continuous one would have had. Anything less and save/resume would quietly
    /// diverge.
    ///
    /// Note this is the *dynamic* pressure — the hydrostatic weight of the
    /// overburden is computed separately and is what melting points respond to.
    pub pressure: &'a mut [f64],
}

/// A mutable view of one grid's state.
pub struct GridView<'a> {
    pub spec: GridSpec,
    /// Micrograms. Mutable now: with a velocity field, matter moves.
    pub mass: &'a mut [i128],
    /// Material index per cell; `MaterialId::VACUUM` for empty.
    pub material: &'a [u16],
    /// Microjoules.
    pub energy: &'a mut [i128],
    /// `None` when the grid holds no fluid — in which case this is Phase 2
    /// exactly as it was, and every Phase 2 test still passes.
    pub fluid: Option<FluidView<'a>>,
}

impl<'a> GridView<'a> {
    pub fn new(
        spec: GridSpec,
        mass: &'a mut [i128],
        material: &'a [u16],
        energy: &'a mut [i128],
    ) -> Self {
        let n = spec.cells();
        assert_eq!(mass.len(), n, "mass column is the wrong length");
        assert_eq!(material.len(), n, "material column is the wrong length");
        assert_eq!(energy.len(), n, "energy column is the wrong length");
        GridView { spec, mass, material, energy, fluid: None }
    }

    pub fn with_fluid(mut self, f: FluidView<'a>) -> Self {
        let n = self.spec.cells();
        assert_eq!(f.mom_x.len(), n);
        assert_eq!(f.mom_y.len(), n);
        assert_eq!(f.mom_z.len(), n);
        assert_eq!(f.pressure.len(), n);
        self.fluid = Some(f);
        self
    }

    #[inline]
    pub fn material_at(&self, i: usize) -> MaterialId {
        MaterialId(self.material[i])
    }

    /// Total energy in the grid. Exact — this is the number the conservation test
    /// watches.
    pub fn total_energy(&self) -> i128 {
        self.energy.iter().sum()
    }

    pub fn total_mass(&self) -> i128 {
        self.mass.iter().sum()
    }

    /// Total momentum, per axis. Exact. Must equal the ledger's impulse.
    pub fn total_momentum(&self) -> [i128; 3] {
        match &self.fluid {
            None => [0; 3],
            Some(f) => [
                f.mom_x.iter().sum(),
                f.mom_y.iter().sum(),
                f.mom_z.iter().sum(),
            ],
        }
    }
}
