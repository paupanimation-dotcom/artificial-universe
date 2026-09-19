//! The planet, validated against physics that was worked out long before software.
//!
//! Phase 2b had Ra_c = 657.5, Phase 3 had van 't Hoff. Phase 4's number-free
//! prediction is the **habitable zone**: the orbital range over which a planet's
//! surface water is liquid rather than ice or vapour. Nothing in this crate draws
//! that band. It falls out of Stefan–Boltzmann setting the temperature and water's
//! phase diagram deciding what that temperature *means* for water. If the engine
//! puts Earth in the liquid-water range, freezes its ocean when we push it out to
//! Mars's orbit, and boils it when we pull it in past Venus's, then the derivation
//! has captured the real physics — and it did so by *using* the exact laws Phase 2
//! already validated, not by adding new ones.

use au_physics::{Material, MaterialRegistry, Phase};
use au_planet::{Planet, PlanetParams, SurfaceState, AU, SOLAR_LUMINOSITY};

const REFERENCE_MATERIALS: &str = include_str!("../../../data/reference_materials.kv");

/// The material registry the physics engine itself uses, loaded from the reference
/// fixture — so a planet's notion of "is water liquid" is the identical phase model
/// that governs a melting ice cube in a Phase 2 cell. An instrument, not content.
fn materials() -> MaterialRegistry {
    let mut m = std::collections::BTreeMap::new();
    for line in REFERENCE_MATERIALS.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some((k, v)) = line.split_once('=') {
            m.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    MaterialRegistry::from_kv(&m).expect("reference materials must parse")
}

fn water(mats: &MaterialRegistry) -> &Material {
    let id = mats.by_name("water").expect("fixture must define water");
    mats.get(id).unwrap()
}

// ═══ GRAVITY: Newton, from mass and radius ══════════════════════════════════

/// Earth's surface gravity comes out at ~9.8 m/s² from its mass and radius alone.
/// Not stored, not declared — computed, and the same g the physics engine uses for
/// pressure and buoyancy.
#[test]
fn earth_surface_gravity_is_about_9_8() {
    let earth = Planet::new(PlanetParams::earth_like());
    let g = earth.params.surface_gravity();
    assert!((g - 9.82).abs() < 0.1, "Earth's gravity should be ~9.8 m/s², got {:.3}", g);
}

/// Gravity scales the right way: a planet of the same density but twice the radius
/// has twice the surface gravity (g ∝ ρ·r). A sanity check that the law, not a
/// table, is doing the work.
#[test]
fn gravity_scales_with_size_at_fixed_density() {
    let base = PlanetParams::earth_like();
    // Double the radius, eight times the mass (same density: M ∝ r³).
    let bigger =
        PlanetParams { radius_m: base.radius_m * 2.0, mass_kg: base.mass_kg * 8.0, ..base };
    let g0 = Planet::new(base).params.surface_gravity();
    let g1 = Planet::new(bigger).params.surface_gravity();
    assert!((g1 / g0 - 2.0).abs() < 1e-9, "doubling radius at fixed density should double g");
}

// ═══ TEMPERATURE: Stefan–Boltzmann, from flux and albedo ════════════════════

/// Earth's equilibrium temperature comes out at ~255 K — the textbook effective
/// temperature, the number every planetary-science course derives. (The warmer
/// *surface* is the greenhouse effect, deliberately not modelled here; the
/// atmosphere at this phase is transparent to its own thermal radiation.) This is
/// Stefan–Boltzmann used globally, the law Phase 2 validated locally against
/// Stefan (1879).
#[test]
fn earth_equilibrium_temperature_is_about_255_k() {
    let earth = Planet::new(PlanetParams::earth_like());
    let t = earth.surface_temperature();
    assert!(
        (t - 255.0).abs() < 5.0,
        "Earth's effective temperature should be ~255 K, got {:.1}",
        t
    );
}

/// The inverse-square law bites: a planet twice as far from its star is colder by
/// a factor of √2 in temperature (T ∝ 1/√d, since flux ∝ 1/d² and T ∝ flux^¼).
#[test]
fn doubling_orbital_distance_cools_by_root_two() {
    let base = PlanetParams::earth_like();
    let far = PlanetParams { orbit_m: base.orbit_m * 2.0, ..base };
    let t0 = Planet::new(base).surface_temperature();
    let t1 = Planet::new(far).surface_temperature();
    let ratio = t0 / t1;
    assert!(
        (ratio - 2f64.sqrt()).abs() < 0.01,
        "twice as far should be √2 cooler, got ratio {:.4}",
        ratio
    );
}

/// Temperature is independent of planet size — the r² cancels in the energy
/// balance. Two planets at the same orbit and albedo have the same mean
/// temperature whether one is a moon and the other a giant.
#[test]
fn temperature_does_not_depend_on_planet_size() {
    let base = PlanetParams::earth_like();
    let tiny =
        PlanetParams { radius_m: base.radius_m * 0.1, mass_kg: base.mass_kg * 0.001, ..base };
    let t0 = Planet::new(base).surface_temperature();
    let t1 = Planet::new(tiny).surface_temperature();
    assert!((t0 - t1).abs() < 1e-6, "temperature should not depend on size");
}

// ═══ THE HABITABLE ZONE: water's phase diagram decides ══════════════════════

/// **The Phase 4 headline.** With Earth's parameters and a surface temperature in
/// the habitable range, water's own phase model — the identical `melt_at`/`boil_at`
/// that governs a Phase 2 ice cube — reports that surface water is **liquid**.
/// Nobody placed an ocean. The star's warmth put the surface in a temperature
/// range, and water answered that in that range it is a liquid.
///
/// An honesty note on the model's current limits: Earth's *effective* temperature
/// is 255 K, below water's 273 K freezing point — without a greenhouse effect a
/// bare-rock Earth freezes, which is a real and famous result (the snowball /
/// faint-young-Sun tension). So this test lowers the albedo until the effective
/// temperature reaches ~288 K, standing in for the greenhouse warming a real
/// atmosphere provides, and isolating the claim under test: *given* a surface
/// temperature, the phase is derived, never painted. The greenhouse itself is
/// later-phase physics (the atmosphere is currently transparent to its own
/// thermal radiation).
#[test]
fn at_earths_orbit_water_is_liquid() {
    let mats = materials();
    let mut p = PlanetParams::earth_like();
    // Solve albedo for T_eff ≈ 288 K: T ∝ (1−A)^¼ ⇒ (1−A) scales by (288/T₀)⁴.
    let want = 288.0_f64;
    let base_t = Planet::new(p).surface_temperature();
    let factor = (want / base_t).powi(4);
    p.albedo = 1.0 - (1.0 - p.albedo) * factor;
    let planet = Planet::new(p);

    let t = planet.surface_temperature();
    assert!((t - 288.0).abs() < 3.0, "setup should put the surface near 288 K, got {:.1}", t);

    let state = planet.surface_state_of(water(&mats));
    assert_eq!(
        state,
        SurfaceState::Liquid,
        "at a habitable surface temperature ({:.0} K), water should be liquid — an ocean, derived not painted",
        t
    );
    assert_eq!(state.phase(), Phase::Liquid);
}

/// Push the planet out toward Mars's orbit and the ocean **freezes** — not because
/// we swapped the ocean for an ice sheet, but because at the colder temperature the
/// star now delivers, water's phase model reports solid. The outer edge of the
/// habitable zone, emergent.
#[test]
fn far_from_the_star_the_ocean_freezes() {
    let mats = materials();
    let p = PlanetParams { orbit_m: 1.6 * AU, ..PlanetParams::earth_like() }; // ~Mars distance
    let planet = Planet::new(p);
    let state = planet.surface_state_of(water(&mats));
    assert_eq!(
        state,
        SurfaceState::Frozen,
        "at 1.6 AU the surface is {:.0} K; water should freeze",
        planet.surface_temperature()
    );
}

/// Pull the planet close enough and the ocean **boils away** to vapour — decided
/// by water's boiling point at the surface pressure, not by fiat. The inner edge
/// of the habitable zone.
///
/// A detail the engine got right and my first version of this test got wrong: at
/// Venus's orbit (0.72 AU) — even at 0.55 AU — the *effective* temperature is only
/// ~345 K, below water's 373 K boiling point at 1 atm. Without a greenhouse, a
/// planet at Venus's distance would have a hot but liquid ocean. Venus is hell
/// because of its runaway greenhouse atmosphere, not its orbit — and this model,
/// whose atmosphere is transparent to thermal radiation, honestly reports that.
/// So the test orbit sits at 0.40 AU, where the bare radiative balance alone
/// (~403 K) exceeds boiling.
#[test]
fn close_to_the_star_the_ocean_vaporises() {
    let mats = materials();
    let p = PlanetParams { orbit_m: 0.40 * AU, ..PlanetParams::earth_like() };
    let planet = Planet::new(p);
    let state = planet.surface_state_of(water(&mats));
    assert_eq!(
        state,
        SurfaceState::Vapor,
        "at 0.40 AU the surface is {:.0} K; water should boil",
        planet.surface_temperature()
    );
}

/// The habitable zone is a *contiguous band*: sweeping orbital distance outward,
/// the phases must go vapour → liquid → frozen monotonically, each a single run —
/// hot inside, cold outside, ocean between. This is the band nobody drew.
///
/// The sweep uses the low-albedo (greenhouse-proxy) planet from the headline test,
/// so a liquid band exists to find; with the bare-rock 255 K Earth the band would
/// sit closer to the star, which is itself physically correct.
#[test]
fn the_habitable_zone_is_a_single_contiguous_band() {
    let mats = materials();
    let w = water(&mats);

    let mut base = PlanetParams::earth_like();
    let base_t = Planet::new(base).surface_temperature();
    base.albedo = 1.0 - (1.0 - base.albedo) * (288.0f64 / base_t).powi(4);

    let mut states = Vec::new();
    let mut d = 0.3;
    while d <= 3.0 {
        let p = PlanetParams { orbit_m: d * AU, ..base };
        states.push((d, Planet::new(p).surface_state_of(w)));
        d += 0.05;
    }

    let order = |s: SurfaceState| match s {
        SurfaceState::Vapor => 0,
        SurfaceState::Liquid => 1,
        SurfaceState::Frozen => 2,
    };
    let mut rank = 0;
    for (d, s) in &states {
        let r = order(*s);
        assert!(
            r >= rank,
            "phases must go vapour→liquid→frozen monotonically outward, but at {:.2} AU it went backwards to {:?}",
            d,
            s
        );
        rank = r;
    }

    assert!(
        states.iter().any(|(_, s)| *s == SurfaceState::Liquid),
        "there should be some orbit where water is liquid"
    );
    assert!(
        states.iter().any(|(_, s)| *s == SurfaceState::Vapor),
        "close in it should be vapour"
    );
    assert!(
        states.iter().any(|(_, s)| *s == SurfaceState::Frozen),
        "far out it should be frozen"
    );
}

/// A brighter star pushes the whole habitable zone outward: flux ∝ L/d², so the
/// same temperature sits at d ∝ √L. Four times the luminosity moves the inner edge
/// out by ~2×. (This is why we look for Earths farther from bright stars.)
#[test]
fn a_brighter_star_pushes_the_habitable_zone_outward() {
    let mats = materials();
    let w = water(&mats);

    let inner_edge = |lum: f64| -> f64 {
        let mut d = 0.2;
        while d <= 6.0 {
            let p = PlanetParams {
                orbit_m: d * AU,
                star_luminosity_w: lum,
                ..PlanetParams::earth_like()
            };
            if Planet::new(p).surface_state_of(w) != SurfaceState::Vapor {
                return d;
            }
            d += 0.01;
        }
        f64::INFINITY
    };

    let dim = inner_edge(SOLAR_LUMINOSITY);
    let bright = inner_edge(4.0 * SOLAR_LUMINOSITY);
    assert!(
        bright > dim * 1.5,
        "a 4× brighter star should push the inner edge out ~2×: {:.2} → {:.2} AU",
        dim,
        bright
    );
}
