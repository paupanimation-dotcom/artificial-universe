//! Terrain relief — the one place the seed shapes structure directly, and the
//! documentation of exactly how far that licence extends.
//!
//! # What this is, honestly
//!
//! The planet's *gross* state is derived by physics (see `planet.rs`): how warm,
//! how heavy, whether water stands liquid. But a real world is not a billiard
//! ball, and the causally honest source of relief — mantle convection pushing up
//! mountains, water carving valleys over geological time — is exactly the
//! deep-time computation Phase 2 proved unreachable (you cannot fast-forward
//! diffusion, and orogeny is diffusion's slowest cousin).
//!
//! So relief is treated as **initial-condition detail**, the same epistemic
//! category as the seed itself: a declared amplitude and wavelength, with the
//! *particular* arrangement of highs and lows drawn deterministically from the
//! world seed. The noise decides only *where the rock stands a little higher*.
//! Everything about what that rock **is** — its temperature, its phase, whether
//! the hollow beside it fills with ocean or ice — is decided by physics, cell by
//! cell, from the planet's derived envelope. A hot world does not get glaciers
//! because the noise felt like it; a hollow below sea level floods because water
//! is liquid *here*, at *this* temperature, and flows to the bottom.
//!
//! When tectonics arrives (it needs the implicit viscous solver, the same one
//! circulation needs), relief becomes an *output* of simulation and this module
//! retires to seeding only the primordial surface those processes then rework.
//! That is the honest trajectory: declared detail today, caused detail tomorrow,
//! and the boundary between them written down rather than blurred.
//!
//! # Determinism and continuity
//!
//! Elevation is a pure function of `(seed, x, y)` — value noise on integer
//! lattices whose corner values come from `Rng::derive(seed, PLANET, key(x), key(y))`,
//! the project's splittable RNG discipline (`f(seed, domain, where)`, nothing
//! global, nothing sequential). Because corners are keyed by *position*, two
//! chunks that share a border sample identical lattice values and the surface is
//! seamless without any communication between them — the same property that lets
//! an untouched world be a function of its seed.
//!
//! Two octaves, smoothstep-interpolated: one at the declared wavelength for the
//! broad shape, one at a quarter wavelength and a quarter amplitude for texture.
//! Nothing here is physical; it is not pretending to be. It is the shape of the
//! dice, declared as such.

use au_core::rng::{Domain, Rng};

/// A deterministic elevation field. Cheap to copy; carries no state.
#[derive(Clone, Copy, Debug)]
pub struct Terrain {
    /// The world seed — the same one everything else derives from.
    pub seed: u64,
    /// Peak-to-mean relief, metres. Elevation stays within ±(amplitude × 1.25)
    /// (the second octave adds a quarter more).
    pub amplitude_m: f64,
    /// Horizontal scale of the broad relief, metres — the distance over which the
    /// landscape changes character.
    pub wavelength_m: f64,
}

impl Terrain {
    pub fn new(seed: u64, amplitude_m: f64, wavelength_m: f64) -> Terrain {
        Terrain { seed, amplitude_m, wavelength_m }
    }

    /// The elevation at a point, metres relative to the reference surface
    /// (positive = above sea level). Pure in `(self, x, y)`.
    pub fn elevation_m(&self, x_m: f64, y_m: f64) -> f64 {
        if self.amplitude_m == 0.0 || self.wavelength_m <= 0.0 {
            return 0.0;
        }
        let broad = self.octave(x_m, y_m, self.wavelength_m, 1);
        let fine = self.octave(x_m, y_m, self.wavelength_m * 0.25, 2);
        self.amplitude_m * (broad + 0.25 * fine)
    }

    /// One octave of value noise in [−1, 1]: bilinear interpolation, smoothstep
    /// eased, between lattice-corner values drawn from the position-keyed RNG.
    fn octave(&self, x_m: f64, y_m: f64, cell_m: f64, octave: u64) -> f64 {
        let fx = x_m / cell_m;
        let fy = y_m / cell_m;
        let x0 = fx.floor();
        let y0 = fy.floor();
        let tx = smooth(fx - x0);
        let ty = smooth(fy - y0);
        let (x0, y0) = (x0 as i64, y0 as i64);

        let v00 = self.corner(x0, y0, octave);
        let v10 = self.corner(x0 + 1, y0, octave);
        let v01 = self.corner(x0, y0 + 1, octave);
        let v11 = self.corner(x0 + 1, y0 + 1, octave);

        let top = v00 + (v10 - v00) * tx;
        let bot = v01 + (v11 - v01) * tx;
        top + (bot - top) * ty
    }

    /// The value at a lattice corner, in [−1, 1] — a fresh RNG stream derived
    /// from (seed, PLANET domain, corner position, octave). Position-keyed, so any
    /// chunk asking about this corner gets the same answer; octave-salted, so the
    /// two octaves are independent fields rather than scaled copies.
    fn corner(&self, x: i64, y: i64, octave: u64) -> f64 {
        // Fold the octave into the y key's high bits: coordinates are far below
        // 2⁴⁸ in practice, so the salt cannot collide with a real position.
        let key_b = (y as u64).wrapping_add(octave << 48);
        let mut rng = Rng::derive(self.seed, Domain::PLANET, x as u64, key_b);
        rng.next_f64() * 2.0 - 1.0
    }
}

/// Smoothstep: eases interpolation so the field's slope is continuous at lattice
/// lines — no creases in the landscape at the seams of the dice.
#[inline]
fn smooth(t: f64) -> f64 {
    t * t * (3.0 - 2.0 * t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn elevation_is_deterministic() {
        let a = Terrain::new(42, 1000.0, 5000.0);
        let b = Terrain::new(42, 1000.0, 5000.0);
        for i in 0..100 {
            let (x, y) = (i as f64 * 137.7, i as f64 * -91.3);
            assert_eq!(a.elevation_m(x, y), b.elevation_m(x, y));
        }
    }

    #[test]
    fn different_seeds_give_different_worlds() {
        let a = Terrain::new(1, 1000.0, 5000.0);
        let b = Terrain::new(2, 1000.0, 5000.0);
        let differs = (0..50).any(|i| {
            let (x, y) = (i as f64 * 313.0, i as f64 * 71.0);
            (a.elevation_m(x, y) - b.elevation_m(x, y)).abs() > 1e-9
        });
        assert!(differs, "two seeds should not share a landscape");
    }

    /// The field is continuous — in particular across the lattice lines where
    /// chunks generated independently must agree. Fine steps never jump.
    #[test]
    fn elevation_is_continuous() {
        let t = Terrain::new(7, 1000.0, 5000.0);
        let mut prev = t.elevation_m(-10_000.0, 333.0);
        let step = 10.0; // metres — far finer than any lattice
        let mut x = -10_000.0 + step;
        while x < 10_000.0 {
            let e = t.elevation_m(x, 333.0);
            assert!(
                (e - prev).abs() < 25.0,
                "a {} m step jumped {:.1} m — a crease at x = {}",
                step,
                (e - prev).abs(),
                x
            );
            prev = e;
            x += step;
        }
    }

    #[test]
    fn elevation_is_bounded_by_the_declared_amplitude() {
        let t = Terrain::new(3, 800.0, 4000.0);
        for i in 0..2000 {
            let (x, y) = (i as f64 * 97.3 - 50_000.0, i as f64 * 41.9 - 30_000.0);
            let e = t.elevation_m(x, y);
            assert!(e.abs() <= 800.0 * 1.25 + 1e-9, "elevation {} exceeds the declared bound", e);
        }
    }
}
