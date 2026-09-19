//! A planet, as a handful of declared physical parameters — and everything that
//! *follows* from them by physical law.
//!
//! # The decision this crate embodies
//!
//! There are three ways to make a planet, and two of them are wrong for this
//! project:
//!
//!   1. **Noise terrain** — a Perlin heightmap painted with biomes. Fast, pretty,
//!      and a flat violation of the Golden Rule: the mountains exist because a
//!      noise function put them there, not because anything caused them. A river
//!      drawn rather than carved. Rejected for the same reason hardcoded molecules
//!      were.
//!
//!   2. **Full physical formation** — simulate the protoplanetary disk, accretion,
//!      radioactive heating, differentiation, tectonics, from dust to world.
//!      Physically honest and utterly uncomputable, worse than the chemistry
//!      cliff. And Phase 2 already proved the deep-time accelerator is inert on
//!      real diffusion, so we could not fast-forward it even if we could afford one
//!      step. Rejected as the cliff.
//!
//!   3. **Procedural from physical parameters** — what this crate does. A planet is
//!      *declared* by its physics: mass, radius, orbital distance, its star's
//!      luminosity, rotation, bulk composition. From those, the planet's
//!      **gross state is derived by physical law**, not painted:
//!
//!        * radius and mass give surface **gravity** (Newton),
//!        * orbital distance and stellar luminosity give the **equilibrium
//!          temperature** (Stefan–Boltzmann — the very law Phase 2 validated),
//!        * temperature and composition decide, through water's own **phase
//!          diagram**, whether the surface volatile is ice, ocean, or vapour —
//!          nobody paints an ocean; the star's warmth and water's boiling point
//!          decide there is one,
//!        * gravity and temperature give the **atmospheric scale height** and the
//!          pressure it falls off with (hydrostatics — also already built).
//!
//! This is the exact analogue of the chemistry decision. There, we declared the
//! *vocabulary* (which elements exist) and let the *outcome* (which reactions fire,
//! what equilibria form) emerge from rules. Here, we declare the *parameters* (how
//! massive, how far from the star) and let the *outcome* (how strong is gravity,
//! is there liquid water) emerge from physics. Declaring that a planet is one Earth
//! mass at one AU is the same kind of act as declaring that carbon exists. It is
//! the opposite of drawing a coastline.
//!
//! # What is derived here, and what is not
//!
//! This module turns parameters into a planet's *gross physical envelope*: gravity,
//! surface temperature, surface pressure, and the phase of its surface volatile.
//! It does **not** produce terrain — that is the generator's job (see the crate
//! root), and even there the fine detail (this chunk's elevation) is procedural
//! *within* the envelope this module sets, never in defiance of it. A world too hot
//! for liquid water does not get an ocean because the seed felt like it.
//!
//! And it does not yet simulate climate, weather, or hydrology — those need
//! circulation, which needs the implicit solver Phase 4's later work will bring.
//! What it gives is the *starting condition* a planet's physics runs forward from:
//! a world in radiative and hydrostatic balance, correct at the largest scale,
//! ready for the finer machinery to act on.

use au_physics::{Material, MaterialRegistry, Phase, P_REF, SIGMA};

/// The gravitational constant, m³·kg⁻¹·s⁻². The one new physical constant this
/// crate needs — Phase 2 never had to weigh a planet.
pub const BIG_G: f64 = 6.674_30e-11;

/// The luminosity of the Sun, W. A convenience reference so stars can be described
/// in solar units; not privileged, just familiar.
pub const SOLAR_LUMINOSITY: f64 = 3.828e26;

/// One astronomical unit, m. Likewise a convenience.
pub const AU: f64 = 1.495_978_707e11;

/// The declared physical parameters of a planet. Everything below is derived from
/// these; these are derived from nothing — they are the planetary analogue of the
/// seed, the givens a universe hands the physics.
#[derive(Clone, Copy, Debug)]
pub struct PlanetParams {
    /// Mass, kg.
    pub mass_kg: f64,
    /// Mean radius, m.
    pub radius_m: f64,
    /// Distance from the star, m. (Circular orbit assumed at this phase; an
    /// eccentric orbit is a seasonal-forcing refinement for later.)
    pub orbit_m: f64,
    /// The star's bolometric luminosity, W.
    pub star_luminosity_w: f64,
    /// Bond albedo, 0..1 — the fraction of incident starlight reflected straight
    /// back to space, never absorbed. Earth's is ~0.29. Higher albedo ⇒ colder.
    pub albedo: f64,
    /// Rotation period, s. Not used for temperature here (a fast rotator evens out
    /// day/night, but the equilibrium *mean* is rotation-independent); carried
    /// because it governs Coriolis forcing, which the circulation work will need.
    pub rotation_s: f64,
    /// Surface pressure of the atmosphere, Pa. Declared for now; a later phase can
    /// derive it from outgassing and escape. Governs, through the phase model,
    /// whether water can be liquid at all (no atmosphere ⇒ water sublimates).
    pub surface_pressure_pa: f64,
}

impl PlanetParams {
    /// An Earth-like planet, for reference and for tests. Real figures.
    pub fn earth_like() -> PlanetParams {
        PlanetParams {
            mass_kg: 5.972e24,
            radius_m: 6.371e6,
            orbit_m: AU,
            star_luminosity_w: SOLAR_LUMINOSITY,
            albedo: 0.29,
            rotation_s: 86_400.0,
            surface_pressure_pa: 101_325.0,
        }
    }

    /// Surface gravity, m·s⁻². Newton: g = G·M / r². The number the physics engine
    /// already knows how to use — it drives hydrostatic pressure and buoyancy.
    pub fn surface_gravity(&self) -> f64 {
        BIG_G * self.mass_kg / (self.radius_m * self.radius_m)
    }

    /// The flux of starlight arriving at the top of the atmosphere, W·m⁻². The
    /// inverse-square law: the star's power spread over a sphere at the orbital
    /// radius. Earth's is the familiar ~1361 W·m⁻² solar constant.
    pub fn stellar_flux(&self) -> f64 {
        self.star_luminosity_w / (4.0 * std::f64::consts::PI * self.orbit_m * self.orbit_m)
    }

    /// The planet's **equilibrium temperature**, K — the temperature at which it
    /// radiates away exactly as much energy as it absorbs from its star.
    ///
    /// This is Stefan–Boltzmann, the same law whose *local* form Phase 2 validated
    /// against Stefan (1879). Here it is used globally: a sphere of radius r
    /// intercepts starlight across its disc (area πr²), keeps the fraction
    /// (1 − albedo), and re-emits across its whole surface (area 4πr²) as a black
    /// body. Setting absorbed = emitted and solving for T:
    ///
    /// ```text
    ///   (1 − A)·F·πr²  =  σ·T⁴·4πr²
    ///   T = [ (1 − A)·F / (4σ) ]^¼
    /// ```
    ///
    /// The r² cancels — a planet's mean temperature does not depend on its size,
    /// only on how much flux it catches and how shiny it is. For Earth's numbers
    /// this gives ~255 K, the textbook effective temperature (the ~33 K warmer
    /// *surface* is the greenhouse effect, which this phase does not model — the
    /// atmosphere here is transparent to its own thermal radiation).
    pub fn equilibrium_temperature(&self) -> f64 {
        let absorbed_per_area = (1.0 - self.albedo) * self.stellar_flux() / 4.0;
        (absorbed_per_area / SIGMA).powf(0.25)
    }

    /// Atmospheric scale height, m — the altitude over which pressure falls by a
    /// factor of e. H = kT / (m·g), written here per unit mean molecular mass so a
    /// caller can supply the atmosphere's composition. Governs how thick the air
    /// column is, and (with surface pressure) the pressure a given altitude feels.
    ///
    /// `mean_molecular_mass_kg` is the mass of one average air molecule (e.g. ~4.8e-26
    /// kg for N₂/O₂ air). Returns 0 for an airless or gravity-less world.
    pub fn scale_height(&self, mean_molecular_mass_kg: f64) -> f64 {
        const K_B: f64 = 1.380_649e-23; // Boltzmann constant, J/K
        let g = self.surface_gravity();
        if g <= 0.0 || mean_molecular_mass_kg <= 0.0 {
            return 0.0;
        }
        K_B * self.equilibrium_temperature() / (mean_molecular_mass_kg * g)
    }
}

/// What a volatile substance *is*, at this planet's surface — decided by physics,
/// not declared.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SurfaceState {
    /// Below the substance's freezing point: ice, rock, frozen.
    Frozen,
    /// Between freezing and boiling: liquid. An ocean, if there is enough of it.
    Liquid,
    /// Above boiling (or, at negligible pressure, sublimated): vapour, atmosphere.
    Vapor,
}

/// A planet: its declared parameters, plus a handle to the material physics that
/// turns them into a state. Immutable after construction.
///
/// The material registry is the same one the physics engine uses — a planet's
/// notion of "is water liquid here" is the *identical* phase model that governs a
/// melting ice cube in a Phase 2 cell. There is one physics, and the planet obeys
/// it.
pub struct Planet {
    pub params: PlanetParams,
}

impl Planet {
    pub fn new(params: PlanetParams) -> Planet {
        Planet { params }
    }

    /// The equilibrium temperature — the surface's mean, in radiative balance.
    pub fn surface_temperature(&self) -> f64 {
        self.params.equilibrium_temperature()
    }

    /// **Is the given surface volatile ice, ocean, or vapour here?** — answered by
    /// the substance's own phase diagram at this planet's temperature and pressure.
    ///
    /// This is the heart of the crate's claim. We do not decide there is an ocean.
    /// We ask water, at the temperature the star sets and the pressure the
    /// atmosphere sets, what state it is in — and water answers, through the exact
    /// `melt_at`/`boil_at` model that governs every Phase 2 cell. Move the planet
    /// out and the answer becomes Frozen; move it in and it becomes Vapor. The
    /// habitable zone is not a band we drew on a map; it is the orbital range over
    /// which this function returns `Liquid`.
    pub fn surface_state_of(&self, volatile: &Material) -> SurfaceState {
        let t = self.surface_temperature();
        let p = self.params.surface_pressure_pa;
        let melt = volatile.melt_at(p);
        let boil = volatile.boil_at(p);
        if t < melt {
            SurfaceState::Frozen
        } else if t > boil {
            SurfaceState::Vapor
        } else {
            SurfaceState::Liquid
        }
    }

    /// Convenience: look a volatile up by name in a registry and report its surface
    /// state. `None` if the material is not declared.
    pub fn surface_state_named(
        &self,
        materials: &MaterialRegistry,
        name: &str,
    ) -> Option<SurfaceState> {
        let id = materials.by_name(name)?;
        Some(self.surface_state_of(materials.get(id).unwrap()))
    }

    /// The hydrostatic pressure at a depth `h` below the surface within a fluid of
    /// the given density: p(h) = p_surface + ρ·g·h. The same relation Phase 2 uses
    /// for a fluid column, applied to a planet's ocean or atmosphere. Provided so
    /// the generator can pressurise deep cells consistently with the physics.
    pub fn pressure_at_depth(&self, depth_m: f64, fluid_density: f64) -> f64 {
        self.params.surface_pressure_pa + fluid_density * self.params.surface_gravity() * depth_m
    }
}

/// The reference pressure the phase model linearises about, re-exported so callers
/// building planet configs can reason about it without reaching into au-physics.
pub const REFERENCE_PRESSURE: f64 = P_REF;

/// The set of phases, re-exported for callers that classify surface states without
/// otherwise depending on au-physics.
pub use au_physics::Phase as MaterialPhase;

impl SurfaceState {
    /// The bulk phase this surface state corresponds to — the bridge back to the
    /// physics `Phase`, so a generator can set a cell's material state to match.
    pub fn phase(self) -> Phase {
        match self {
            SurfaceState::Frozen => Phase::Solid,
            SurfaceState::Liquid => Phase::Liquid,
            SurfaceState::Vapor => Phase::Gas,
        }
    }
}
