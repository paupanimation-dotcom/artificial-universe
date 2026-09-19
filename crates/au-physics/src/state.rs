//! Temperature is not a thing you set. It is a thing you get.
//!
//! This file is small and it is the philosophical heart of Phase 2.
//!
//! PROJECT_VISION's Fundamental Rule is that everything exists because something
//! caused it. Applied here, that has a sharp and slightly uncomfortable
//! consequence: **temperature cannot be a stored variable.** If it were, then
//! "the rock is hot" would be a fact the designer asserted, and heating
//! something would mean *writing a number*. Cause and effect would run
//! backwards.
//!
//! What actually exists is *energy in matter*. Temperature is what emerges when
//! you ask how much. So:
//!
//! ```text
//! state = f(mass, energy, material, pressure, volume)
//! ```
//!
//! You cannot set the temperature. You add energy, and the temperature is
//! whatever the substance decides it is. Sometimes — and this is the good part —
//! it decides *nothing at all*, and melts instead.
//!
//! # Latent heat, and why it matters more than it looks
//!
//! Pour energy into ice at 273 K and the temperature does not move. It sits
//! there, absorbing 334 kJ/kg, while the ice becomes water. Only then does it
//! start warming again.
//!
//! That plateau is the reason Earth has a climate rather than a temperature
//! swing. It is a thermostat with no thermostat in it — a stabilising feedback
//! that *falls out of the thermodynamics* rather than being designed. It is
//! precisely the kind of emergent regulation this whole project is a bet on, and
//! it costs about fifteen lines.

use crate::material::{Material, Phase};

/// Everything derivable about a cell right now.
///
/// Computed fresh every tick from the exact integers. Never stored, never
/// accumulated, so float error has nowhere to breed.
#[derive(Clone, Copy, Debug)]
pub struct CellState {
    pub temperature: f64,
    pub phase: Phase,
    /// 0 = fully the lower phase, 1 = fully the upper phase. Nonzero only while
    /// a transition is in progress — this is the cell mid-melt.
    pub transition: f64,
    pub density: f64,
    pub conductivity: f64,
    /// Heat capacity of the whole cell, J/K. Not specific: total.
    pub heat_capacity: f64,
    /// Thermal diffusivity, m²/s. The number that sets the timestep.
    pub diffusivity: f64,
}

impl CellState {
    /// A cell with nothing in it.
    pub const VACUUM: CellState = CellState {
        temperature: 0.0,
        phase: Phase::Gas,
        transition: 0.0,
        density: 0.0,
        conductivity: 0.0,
        heat_capacity: 0.0,
        diffusivity: 0.0,
    };
}

#[inline]
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

/// Derive everything from what actually exists.
///
/// `mass_kg` and `energy_j` come from the exact integer state. `pressure_pa`
/// comes from gravity acting on the matter above. `volume_m3` is the cell.
///
/// Energy is measured from a zero at 0 K in the solid phase. A negative energy
/// is unphysical and means something upstream leaked — the debug assertion is
/// there because that bug must be found in the tick it happens, not in Phase 9.
pub fn derive(mass_kg: f64, energy_j: f64, mat: &Material, pressure_pa: f64, volume_m3: f64) -> CellState {
    if mass_kg <= 0.0 || volume_m3 <= 0.0 {
        return CellState::VACUUM;
    }
    debug_assert!(
        energy_j >= 0.0,
        "cell has negative energy ({} J) — something is leaking",
        energy_j
    );

    let density = mass_kg / volume_m3;
    let tm = mat.melt_at(pressure_pa);
    let tb = mat.boil_at(pressure_pa);

    // The four energy thresholds. Everything below is just deciding which of the
    // five regimes we are in.
    let e_melt_start = mass_kg * mat.c_solid * tm; //           fully solid, at the melting point
    let e_melt_end = e_melt_start + mass_kg * mat.latent_fusion; // fully liquid, still at the melting point
    let e_boil_start = e_melt_end + mass_kg * mat.c_liquid * (tb - tm);
    let e_boil_end = e_boil_start + mass_kg * mat.latent_vapor;

    let e = energy_j.max(0.0);

    let (temperature, phase, transition, k, c) = if e < e_melt_start {
        // Solid, warming.
        let t = e / (mass_kg * mat.c_solid);
        (t, Phase::Solid, 0.0, mat.k_solid, mat.c_solid)
    } else if e < e_melt_end {
        // MELTING. Energy pours in; temperature does not move. The plateau.
        let f = (e - e_melt_start) / (mass_kg * mat.latent_fusion);
        let phase = if f < 0.5 { Phase::Solid } else { Phase::Liquid };
        (
            tm,
            phase,
            f,
            lerp(mat.k_solid, mat.k_liquid, f),
            lerp(mat.c_solid, mat.c_liquid, f),
        )
    } else if e < e_boil_start {
        // Liquid, warming.
        let t = tm + (e - e_melt_end) / (mass_kg * mat.c_liquid);
        (t, Phase::Liquid, 0.0, mat.k_liquid, mat.c_liquid)
    } else if e < e_boil_end {
        // Boiling. The plateau again, an order of magnitude wider.
        let f = (e - e_boil_start) / (mass_kg * mat.latent_vapor);
        let phase = if f < 0.5 { Phase::Liquid } else { Phase::Gas };
        (
            tb,
            phase,
            f,
            lerp(mat.k_liquid, mat.k_gas, f),
            lerp(mat.c_liquid, mat.c_gas, f),
        )
    } else {
        // Gas, warming.
        let t = tb + (e - e_boil_end) / (mass_kg * mat.c_gas);
        (t, Phase::Gas, 0.0, mat.k_gas, mat.c_gas)
    };

    CellState {
        temperature,
        phase,
        transition,
        density,
        conductivity: k,
        heat_capacity: mass_kg * c,
        diffusivity: if density > 0.0 && c > 0.0 { k / (density * c) } else { 0.0 },
    }
}

/// The inverse, for *setting up* a world — never for running one.
///
/// Tests and (in Phase 4) planet generation need to say "make this cell 300 K".
/// That is a legitimate initial condition. It is not a legitimate *operation*:
/// nothing in the running simulation may call this, because heating something by
/// assigning its temperature is exactly the causal inversion this file exists to
/// prevent. It lives here, named awkwardly, so that its use is conspicuous.
pub fn energy_for_temperature(
    mass_kg: f64,
    target_k: f64,
    mat: &Material,
    pressure_pa: f64,
) -> f64 {
    let tm = mat.melt_at(pressure_pa);
    let tb = mat.boil_at(pressure_pa);
    if target_k <= tm {
        mass_kg * mat.c_solid * target_k
    } else if target_k <= tb {
        mass_kg * (mat.c_solid * tm + mat.latent_fusion + mat.c_liquid * (target_k - tm))
    } else {
        mass_kg
            * (mat.c_solid * tm
                + mat.latent_fusion
                + mat.c_liquid * (tb - tm)
                + mat.latent_vapor
                + mat.c_gas * (target_k - tb))
    }
}
