//! The planetary chunk generator — where a declared planet becomes stored matter.
//!
//! This is the third of the "where X meets storage" files (`physics_columns`,
//! `chem_columns`, now this), and the first *generator* the project has that is
//! not a test fixture. Until now every active chunk was hand-seeded by a test or
//! demo; from Phase 4 on, a world can declare a planet in config and the chunks
//! fill themselves — matter placed as a pure function of `(seed, chunk coord,
//! planet parameters)`, honouring the contract that an untouched world is a
//! function of its seed.
//!
//! # Who decides what, cell by cell
//!
//! For every cell the generator asks two questions, and the *sources* of the two
//! answers are the whole point:
//!
//!   * **Where does the rock stand?** — the terrain field answers: seeded relief,
//!     the honestly-declared initial-condition detail (`au_planet::terrain`).
//!
//!   * **What is everything's state?** — physics answers, and only physics. A cell
//!     below the rock line is rock at the temperature the star sets. A hollow
//!     below sea level fills with the planet's volatile — and whether that
//!     volatile is an ocean or an ice sheet is *not chosen here*. The generator
//!     writes mass and energy; the phase is whatever `au_physics::derive` reports
//!     when anyone looks, at the energy the star's derived temperature implies. A
//!     cold planet's "ocean" cells come out solid — a world of ice — because the
//!     physics says so, not because the generator painted glaciers. Move the same
//!     config closer to the star, change nothing else, and the identical cells
//!     come out liquid.
//!
//! That division — seed decides *where*, physics decides *what* — is the whole
//! licence the Golden Rule grants, and this file stays carefully inside it.
//!
//! # The vertical picture
//!
//! World z is altitude in metres, sea level at z = 0. A cell's centre below the
//! local terrain elevation is rock; between rock and sea level it is volatile;
//! above both it is vacuum (a real atmosphere as stored gas cells needs the
//! circulation work — at this phase the "air" is the declared surface pressure the
//! phase model feels, not resident matter). Cells are filled whole; partial-cell
//! coastlines are a refinement that belongs with rendering, not physics.
//!
//! # Determinism
//!
//! The generator ignores the per-chunk RNG handed to it: relief must agree across
//! chunk borders, so it derives its randomness position-keyed inside the terrain
//! field instead (`Rng::derive(seed, PLANET, x, y)` — same discipline, different
//! keying). Two calls with the same seed and coordinate produce byte-identical
//! chunks; the world-hash tests hold the generator to that.

use au_core::Rng;
use au_data::chunk::{Chunk, ChunkCoord, ChunkGenerator, Lod};
use au_physics::{energy_for_temperature, Energy, GridSpec, Mass, MaterialId, MaterialRegistry};
use au_planet::{Planet, PlanetParams, Terrain, AU, SOLAR_LUMINOSITY};

use crate::config::Config;
use crate::physics_columns::{install, PHYS_ENERGY, PHYS_MASS, PHYS_MATERIAL};

/// A generator that fills chunks with a declared planet's matter.
pub struct PlanetGenerator {
    pub planet: Planet,
    pub terrain: Terrain,
    spec: GridSpec,
    rock: MaterialId,
    ocean: MaterialId,
    rock_material: au_physics::Material,
    ocean_material: au_physics::Material,
    rock_density: f64,
    ocean_density: f64,
    /// The planet's derived surface temperature, cached — one Stefan–Boltzmann
    /// evaluation, used for every cell.
    surface_k: f64,
    n_chem_species: usize,
}

impl PlanetGenerator {
    /// Build from config, or `None` if no planet is declared. Panics on a planet
    /// declared with materials the world does not define — a planet made of
    /// nothing is a config bug, not a condition to limp past.
    pub fn from_config(
        seed: u64,
        c: &Config,
        materials: &MaterialRegistry,
        n_chem_species: usize,
    ) -> Option<PlanetGenerator> {
        if !c.bool_or("planet.enabled", false) {
            return None;
        }
        let params = PlanetParams {
            mass_kg: c.f64_or("planet.mass_kg", 5.972e24),
            radius_m: c.f64_or("planet.radius_m", 6.371e6),
            orbit_m: c.f64_or("planet.orbit_au", 1.0) * AU,
            star_luminosity_w: c.f64_or("planet.star_luminosity_solar", 1.0) * SOLAR_LUMINOSITY,
            albedo: c.f64_or("planet.albedo", 0.29),
            rotation_s: c.f64_or("planet.rotation_s", 86_400.0),
            surface_pressure_pa: c.f64_or("planet.surface_pressure_pa", 101_325.0),
        };
        let planet = Planet::new(params);
        let terrain = Terrain::new(
            seed,
            c.f64_or("planet.relief_m", 60.0),
            c.f64_or("planet.relief_wavelength_m", 400.0),
        );
        let spec = GridSpec::new(
            c.u64_or("physics.grid.nx", 1) as usize,
            c.u64_or("physics.grid.ny", 1) as usize,
            c.u64_or("physics.grid.nz", 1) as usize,
            c.f64_or("physics.grid.cell_m", 10.0),
        );
        let rock_name = c.str_or("planet.rock", "silicate");
        let ocean_name = c.str_or("planet.ocean", "water");
        let rock = materials
            .by_name(&rock_name)
            .unwrap_or_else(|| panic!("planet.rock '{}' is not a declared material", rock_name));
        let ocean = materials
            .by_name(&ocean_name)
            .unwrap_or_else(|| panic!("planet.ocean '{}' is not a declared material", ocean_name));
        let surface_k = planet.surface_temperature();

        let rock_material = materials.get(rock).unwrap().clone();
        let ocean_material = materials.get(ocean).unwrap().clone();
        Some(PlanetGenerator {
            planet,
            terrain,
            spec,
            rock,
            ocean,
            rock_material,
            ocean_material,
            rock_density: c.f64_or("planet.rock_density", 3000.0),
            ocean_density: c.f64_or("planet.ocean_density", 1000.0),
            surface_k,
            n_chem_species,
        })
    }

    /// The materials this generator will place. Exposed for tests and demos that
    /// want to interrogate what the physics reports about the generated cells.
    pub fn rock_id(&self) -> MaterialId {
        self.rock
    }
    pub fn ocean_id(&self) -> MaterialId {
        self.ocean
    }
}

impl ChunkGenerator for PlanetGenerator {
    fn generate(&self, coord: ChunkCoord, _rng: &mut Rng) -> Chunk {
        let mut chunk = Chunk::new(coord, Lod::Full);
        let cells = self.spec.cells();
        install(&mut chunk.columns, cells);
        if self.n_chem_species > 0 {
            crate::chem_columns::install(&mut chunk.columns, cells, self.n_chem_species);
        }

        let cell_m = self.spec.cell_m;
        let vol = self.spec.cell_volume();
        // The chunk's world origin, metres: chunk coordinate × chunk extent.
        let ox = coord.x as f64 * self.spec.nx as f64 * cell_m;
        let oy = coord.y as f64 * self.spec.ny as f64 * cell_m;
        let oz = coord.z as f64 * self.spec.nz as f64 * cell_m;

        let p_surface = self.planet.params.surface_pressure_pa;
        // Precompute per-material fill values once — every rock cell is identical,
        // every volatile cell is identical; only the *arrangement* varies.
        let rock_kg = self.rock_density * vol;
        let ocean_kg = self.ocean_density * vol;

        // Fill.
        let mut mass = vec![0i128; cells];
        let mut energy = vec![0i128; cells];
        let mut material = vec![MaterialId::VACUUM.0; cells];

        for z in 0..self.spec.nz {
            let z_m = oz + (z as f64 + 0.5) * cell_m; // cell-centre altitude
            for y in 0..self.spec.ny {
                let y_m = oy + (y as f64 + 0.5) * cell_m;
                for x in 0..self.spec.nx {
                    let x_m = ox + (x as f64 + 0.5) * cell_m;
                    let i = self.spec.idx(x, y, z);

                    let elev = self.terrain.elevation_m(x_m, y_m);
                    if z_m < elev {
                        // Rock, at the temperature the star sets.
                        material[i] = self.rock.0;
                        mass[i] = Mass::from_kg(rock_kg).0;
                        energy[i] = Energy::from_joules(self.rock_energy(rock_kg, p_surface)).0;
                    } else if z_m < 0.0 {
                        // Below sea level: the volatile. Whether these cells ARE an
                        // ocean or an ice sheet is not decided here — the energy
                        // written implies the star's temperature, and the phase is
                        // whatever the physics derives from it.
                        material[i] = self.ocean.0;
                        mass[i] = Mass::from_kg(ocean_kg).0;
                        energy[i] = Energy::from_joules(self.ocean_energy(ocean_kg, p_surface)).0;
                    }
                    // else: vacuum — already zeroed.
                }
            }
        }

        chunk.columns.get_mut::<i128>(PHYS_MASS).unwrap().as_mut_slice().copy_from_slice(&mass);
        chunk.columns.get_mut::<i128>(PHYS_ENERGY).unwrap().as_mut_slice().copy_from_slice(&energy);
        chunk
            .columns
            .get_mut::<u16>(PHYS_MATERIAL)
            .unwrap()
            .as_mut_slice()
            .copy_from_slice(&material);
        chunk
    }
}

impl PlanetGenerator {
    // The generator owns copies of its two materials (they are small value
    // types), so cell fill needs no registry access — the world keeps sole
    // ownership of the registry, and generation stays a pure function.
    fn rock_energy(&self, kg: f64, pressure: f64) -> f64 {
        energy_for_temperature(kg, self.surface_k, &self.rock_material, pressure)
    }
    fn ocean_energy(&self, kg: f64, pressure: f64) -> f64 {
        energy_for_temperature(kg, self.surface_k, &self.ocean_material, pressure)
    }
}
