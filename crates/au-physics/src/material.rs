//! Bulk materials — and the line between Phase 2 and Phase 3.
//!
//! ARCHITECTURE gives Physics a "Materials" section (solid/liquid/gas, density,
//! temperature effects) *and* gives Matter its own layer (hydrogen, carbon,
//! oxygen...). Those are different things and conflating them would collapse two
//! layers into one:
//!
//! * **Phase 2 (here):** what a substance *does*. It has a heat capacity, it
//!   conducts, it melts, it boils. These are the laws.
//! * **Phase 3 (Matter):** what a substance *is*. Elements, bonds, molecules.
//!   These are the causes.
//!
//! So this file ships the *model* and no substances. The registry starts empty
//! and stays empty, because in Phase 2 there is no such thing as iron — there is
//! only "a thing with these bulk properties". Phase 3's job will be to *derive*
//! entries in this registry from element composition, and at that point a
//! material stops being a hand-written row and starts being a consequence.
//!
//! `data/reference_materials.kv` exists, but it is **not content and the engine
//! never loads it.** It is a validation fixture: real measured Earth substances,
//! used by the tests to check that this engine reproduces reality (a latent-heat
//! plateau at the right temperature, a radiative equilibrium at the right
//! temperature). PROJECT_VISION permits exactly this — "Earth is only a
//! scientific reference for how these systems work."
//!
//! # Note on what is *not* here
//!
//! There is no `density` field. Density is `mass / volume` — it is *derived*,
//! and a material that carried its own density would be asserting a fact the
//! simulation is supposed to compute. If a cell is denser than its substance
//! "should" be, that is compression, and compression is a physical result, not a
//! constant.

use std::collections::BTreeMap;

/// Reference pressure for the material table: 1 standard atmosphere.
/// Melting and boiling points are quoted *at* this pressure and shift away from
/// it via the Clausius–Clapeyron slopes below.
pub const P_REF: f64 = 101_325.0;

/// Stefan–Boltzmann constant, W·m⁻²·K⁻⁴.
pub const SIGMA: f64 = 5.670_374_419e-8;

/// State of matter.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
#[repr(u8)]
pub enum Phase {
    Solid = 0,
    Liquid = 1,
    Gas = 2,
}

impl Phase {
    pub fn name(self) -> &'static str {
        match self {
            Phase::Solid => "solid",
            Phase::Liquid => "liquid",
            Phase::Gas => "gas",
        }
    }
}

/// Index into the [`MaterialRegistry`]. Stored per cell as a `u16`.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub struct MaterialId(pub u16);

impl MaterialId {
    /// Nothing. A cell with no material is vacuum: it holds no mass, conducts
    /// nothing, and radiates nothing.
    pub const VACUUM: MaterialId = MaterialId(u16::MAX);
}

/// The bulk thermodynamic behaviour of a substance.
#[derive(Clone, Debug)]
pub struct Material {
    pub name: String,

    /// Specific heat capacity, J·kg⁻¹·K⁻¹, per phase.
    pub c_solid: f64,
    pub c_liquid: f64,
    pub c_gas: f64,

    /// Thermal conductivity, W·m⁻¹·K⁻¹, per phase.
    pub k_solid: f64,
    pub k_liquid: f64,
    pub k_gas: f64,

    /// Transition temperatures at [`P_REF`], in kelvin.
    pub melt_k: f64,
    pub boil_k: f64,

    /// Latent heats, J·kg⁻¹. This is where the interesting physics lives:
    /// energy that goes into *changing state* rather than raising temperature.
    /// It is why oceans are a thermostat, and why a planet with a volatile can
    /// hold a stable climate at all.
    pub latent_fusion: f64,
    pub latent_vapor: f64,

    /// Clausius–Clapeyron slopes, K·Pa⁻¹ — how the transition temperatures move
    /// with pressure. Linearised, which is a lie at extreme pressures and a good
    /// one at planetary ones.
    ///
    /// Not a decoration: this is why Earth's inner core is *solid* despite being
    /// hotter than the liquid outer core. Structure that arises from pressure,
    /// not from anyone deciding a planet should have a core.
    pub dtm_dp: f64,
    pub dtb_dp: f64,

    /// Thermal emissivity, 0..1. How well an exposed surface radiates.
    pub emissivity: f64,

    /// Dynamic viscosity, Pa·s, per fluid phase. A solid does not flow at all —
    /// it is treated as rigid, not as a very thick liquid, because "very thick"
    /// would mean a viscous timestep of ~10⁻¹³ s and the simulation would stop.
    pub mu_liquid: f64,
    pub mu_gas: f64,

    /// Volumetric thermal expansion coefficient, K⁻¹.
    ///
    /// **The single number that makes convection possible.** Without it, heating
    /// a parcel changes its temperature and nothing else — it is exactly as heavy
    /// as its neighbours, so it does not rise, so nothing ever stirs. With it, a
    /// hot parcel is lighter than the cold fluid around it, and a hot parcel that
    /// is lighter than its surroundings is a parcel that goes *up*.
    ///
    /// Everything downstream — mantle overturn, ocean circulation, weather, the
    /// hydrothermal vent that may one day matter more than any of them — is
    /// downstream of this coefficient being nonzero.
    pub alpha: f64,
}

impl Material {
    #[inline]
    pub fn specific_heat(&self, phase: Phase) -> f64 {
        match phase {
            Phase::Solid => self.c_solid,
            Phase::Liquid => self.c_liquid,
            Phase::Gas => self.c_gas,
        }
    }

    /// Viscosity by phase. A solid returns `None` — it is rigid, not viscous.
    #[inline]
    pub fn viscosity(&self, phase: Phase) -> Option<f64> {
        match phase {
            Phase::Solid => None,
            Phase::Liquid => Some(self.mu_liquid),
            Phase::Gas => Some(self.mu_gas),
        }
    }

    #[inline]
    pub fn conductivity(&self, phase: Phase) -> f64 {
        match phase {
            Phase::Solid => self.k_solid,
            Phase::Liquid => self.k_liquid,
            Phase::Gas => self.k_gas,
        }
    }

    /// Melting point at a given pressure.
    #[inline]
    pub fn melt_at(&self, pressure: f64) -> f64 {
        self.melt_k + self.dtm_dp * (pressure - P_REF)
    }

    /// Boiling point at a given pressure. Clamped never to fall below melting:
    /// beyond the triple point the linearisation stops making sense, and a
    /// substance that "boils below its melting point" would send the state
    /// derivation into nonsense.
    #[inline]
    pub fn boil_at(&self, pressure: f64) -> f64 {
        (self.boil_k + self.dtb_dp * (pressure - P_REF)).max(self.melt_at(pressure))
    }

    /// Thermal diffusivity, m²·s⁻¹ — how fast heat *spreads*, as opposed to how
    /// fast it flows. This is the number that decides the timestep, and
    /// therefore the number that decides whether deep time is reachable.
    #[inline]
    pub fn diffusivity(&self, phase: Phase, density: f64) -> f64 {
        if density <= 0.0 {
            return 0.0;
        }
        let c = self.specific_heat(phase);
        if c <= 0.0 {
            return 0.0;
        }
        self.conductivity(phase) / (density * c)
    }
}

/// The substances that exist.
///
/// Empty by default and empty on ship. Phase 3 will fill it — not by hand, but
/// by deriving bulk properties from molecular structure. Until then it is filled
/// only by tests, which is exactly the right amount of content for a phase whose
/// job is to define laws.
#[derive(Clone, Debug, Default)]
pub struct MaterialRegistry {
    materials: Vec<Material>,
}

impl MaterialRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, m: Material) -> MaterialId {
        assert!(self.materials.len() < u16::MAX as usize, "too many materials");
        self.materials.push(m);
        MaterialId((self.materials.len() - 1) as u16)
    }

    #[inline]
    pub fn get(&self, id: MaterialId) -> Option<&Material> {
        self.materials.get(id.0 as usize)
    }

    pub fn by_name(&self, name: &str) -> Option<MaterialId> {
        self.materials
            .iter()
            .position(|m| m.name == name)
            .map(|i| MaterialId(i as u16))
    }

    pub fn len(&self) -> usize {
        self.materials.len()
    }
    pub fn is_empty(&self) -> bool {
        self.materials.is_empty()
    }
    pub fn iter(&self) -> impl Iterator<Item = (MaterialId, &Material)> {
        self.materials
            .iter()
            .enumerate()
            .map(|(i, m)| (MaterialId(i as u16), m))
    }

    /// Parse from `key = value` pairs.
    ///
    /// Config-driven rather than code-driven so that a saved world reloads with
    /// exactly the substances it was computed with. A material table that
    /// differed between save and load would silently change the physics of a
    /// world mid-history.
    pub fn from_kv(kv: &BTreeMap<String, String>) -> Result<MaterialRegistry, String> {
        let count: usize = kv
            .get("material.count")
            .map(|s| s.parse().map_err(|_| "material.count is not an integer".to_string()))
            .transpose()?
            .unwrap_or(0);

        let mut reg = MaterialRegistry::new();
        for i in 0..count {
            let g = |k: &str| -> Result<f64, String> {
                let key = format!("material.{}.{}", i, k);
                kv.get(&key)
                    .ok_or_else(|| format!("missing {}", key))?
                    .parse()
                    .map_err(|_| format!("{} is not a number", key))
            };
            let name = kv
                .get(&format!("material.{}.name", i))
                .cloned()
                .unwrap_or_else(|| format!("material_{}", i));

            reg.add(Material {
                name,
                c_solid: g("c_solid")?,
                c_liquid: g("c_liquid")?,
                c_gas: g("c_gas")?,
                k_solid: g("k_solid")?,
                k_liquid: g("k_liquid")?,
                k_gas: g("k_gas")?,
                melt_k: g("melt_k")?,
                boil_k: g("boil_k")?,
                latent_fusion: g("latent_fusion")?,
                latent_vapor: g("latent_vapor")?,
                dtm_dp: g("dtm_dp")?,
                dtb_dp: g("dtb_dp")?,
                emissivity: g("emissivity")?,
                mu_liquid: g("mu_liquid").unwrap_or(0.0),
                mu_gas: g("mu_gas").unwrap_or(0.0),
                alpha: g("alpha").unwrap_or(0.0),
            });
        }
        Ok(reg)
    }
}
