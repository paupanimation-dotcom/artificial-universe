//! Vesicles — when a bag has enough surface to become two bags.
//!
//! # What is missing, precisely
//!
//! Phase 5f gave the world compartments and Phase 5h gave compartments a reason
//! to exist: in a flow, a barrier keeps its contents, by a factor of 457 in the
//! test that measures it. But a compartment that merely *persists* is a rock. It
//! has no descendants, so nothing about it can be inherited, so the difference
//! between a good barrier and a bad one cannot accumulate. Selection needs
//! something that makes more of itself and passes on what it was.
//!
//! This module is the first half of that: **when does a bag divide, and what do
//! its daughters get?** It decides neither where bags are nor how many there are
//! — that is storage, and it comes next. Everything here is a pure function of a
//! composition and its surroundings.
//!
//! # Nobody chooses a size threshold
//!
//! The tempting implementation is a rule: divide when the contents exceed N.
//! That is a designer's number, it would be the most load-bearing constant in
//! the engine, and the Golden Rule forbids exactly this kind of thing. So the
//! criterion is taken from geometry, where it has been sitting all along.
//!
//! A closed bag has two independent quantities:
//!
//!   * an **area**, set by how many amphiphiles it has assembled —
//!     `A = N · area_per_molecule`, the same relation Phase 5f already used to
//!     ask whether a cell could close a bag at all;
//!   * a **volume**, set by osmosis. A membrane passes water and not solute, so
//!     the bag swells or shrinks until the concentration inside matches the
//!     concentration outside: `V = n_internal / c_external`. Not a knob. A
//!     consequence of what the bag contains and what surrounds it.
//!
//! Those two are made by different chemistry and they do not have to keep pace.
//! Their ratio is the **reduced volume**
//!
//! ```text
//!     v = V / V_sphere(A) = 3·V·√(4π) / A^(3/2)
//! ```
//!
//! which is 1 for a sphere and less than 1 for anything floppier, and it is the
//! standard order parameter of vesicle shape (Seifert, Berndl & Lipowsky, 1991).
//! Two facts about it decide everything here, and both are arithmetic:
//!
//! **v > 1 is impossible.** The sphere encloses the most volume per unit area
//! there is, so a bag whose osmotic volume exceeds its spherical maximum cannot
//! hold itself together. It bursts. Real bilayers tolerate a few percent of
//! areal strain first — about 3%, which is where `lysis_strain` comes from —
//! and then they fail.
//!
//! **v ≤ 1/√2 admits two spheres.** Take area `A` and ask when it is enough to
//! close *two* equal spheres holding the same total volume. Two spheres of
//! radius r have area `8πr²` and volume `(8/3)πr³`, and substituting gives
//! exactly `v = 1/√2 ≈ 0.7071`. Above it, symmetric fission would need more
//! membrane than the bag owns; below it, the bag has surface to spare and two
//! daughters are geometrically available.
//!
//! `1/√2` is not a tuned parameter. It is what `√π/√(2π)` equals.
//!
//! # The trade-off this creates, which nobody wrote down
//!
//! A protocell now sits between two failures. Make osmolytes faster than
//! membrane and `v` climbs toward 1 and it bursts. Make membrane faster than
//! osmolytes and `v` falls to 1/√2 and it divides. **Whether a bag divides,
//! bursts or persists is decided by the ratio of two rates in its own
//! chemistry** — not by its size, not by a timer, and not by anything in this
//! file. That is the shape of answer this project exists to produce.
//!
//! # A caveat worth stating plainly
//!
//! `1/√2` is the *geometric* bound: the point at which two daughters become
//! possible, not the point at which they become energetically favourable. The
//! real budding transition is set by bending rigidity and spontaneous curvature
//! and sits at a higher reduced volume, varying with the lipid. What is
//! implemented here is the hard bound, which is exact and needs no material
//! parameters; the energetics would move the line, not remove it.

use au_core::Rng;

/// One-over-root-two: the reduced volume at which a vesicle's area is exactly
/// enough to close two equal spheres of the same total volume.
///
/// Derived, not chosen. See the module docs.
pub const FISSION_REDUCED_VOLUME: f64 = std::f64::consts::FRAC_1_SQRT_2;

/// What geometry permits a bag to do.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Fate {
    /// Osmotic volume exceeds what the membrane can enclose, past the strain a
    /// bilayer tolerates. The bag fails.
    Lysis,
    /// Enough area for one sphere and not two. It persists.
    Intact,
    /// Enough area for two equal spheres. Division is geometrically available.
    Fission,
}

/// The shape parameters of a bag, all derived from what it holds.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Shape {
    /// Membrane area, m².
    pub area_m2: f64,
    /// Osmotic volume, m³.
    pub volume_m3: f64,
    /// Reduced volume: 1 for a sphere, less for anything floppier.
    pub reduced_volume: f64,
}

/// Geometry of a bag holding `membrane` assembled amphiphiles and `osmolytes`
/// free internal solute molecules, surrounded by a medium at `c_external`
/// molecules·m⁻³.
///
/// `area_per_molecule_m2` is the same headgroup area Phase 5f uses to ask
/// whether a cell can close a bag at all — one constant, one meaning.
pub fn shape(
    membrane: i128,
    osmolytes: i128,
    c_external: f64,
    area_per_molecule_m2: f64,
) -> Option<Shape> {
    if membrane <= 0 || c_external <= 0.0 || area_per_molecule_m2 <= 0.0 {
        return None;
    }
    let area = membrane as f64 * area_per_molecule_m2;
    // Osmotic balance: the membrane passes water, not solute, so the bag swells
    // until what is inside is as dilute as what is outside.
    let volume = (osmolytes.max(0) as f64) / c_external;
    let v = reduced_volume(area, volume)?;
    Some(Shape { area_m2: area, volume_m3: volume, reduced_volume: v })
}

/// `v = 3V√(4π) / A^(3/2)` — 1 for a sphere, less for anything else.
pub fn reduced_volume(area_m2: f64, volume_m3: f64) -> Option<f64> {
    if !(area_m2 > 0.0) || volume_m3 < 0.0 || !area_m2.is_finite() || !volume_m3.is_finite() {
        return None;
    }
    Some(3.0 * volume_m3 * (4.0 * std::f64::consts::PI).sqrt() / area_m2.powf(1.5))
}

/// What geometry permits, given a bag's shape.
///
/// `lysis_strain` is the fractional areal stretch a bilayer survives before it
/// fails — about 0.03 for real lipids. Expressed as strain rather than as a
/// reduced-volume cutoff because that is the quantity experiments measure.
pub fn fate(v: f64, lysis_strain: f64) -> Fate {
    // A bag stretched by areal strain ε encloses (1+ε)^(3/2) times the volume,
    // so tolerating strain ε means tolerating v up to that factor.
    let ceiling = (1.0 + lysis_strain.max(0.0)).powf(1.5);
    if v > ceiling {
        Fate::Lysis
    } else if v <= FISSION_REDUCED_VOLUME {
        Fate::Fission
    } else {
        Fate::Intact
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Partition — what the daughters get
// ═══════════════════════════════════════════════════════════════════════════

/// Split a bag's contents between two daughters.
///
/// # Why this is the engine's first random physics, and why that is not a
/// betrayal of determinism
///
/// Every physical process in this engine so far has had a mean-field answer.
/// Partition does not. When a bag pinches in two, which side a given molecule
/// ends up on is decided by where it happened to be, and there is no
/// deterministic function of the parent that gives the answer — the honest
/// model is a coin per molecule, weighted by the daughters' volumes.
///
/// Phase 5b met this question with diffusion below the quantisation floor and
/// deferred it, on the grounds that introducing randomness was a real decision
/// about the determinism story rather than a rounding detail. This is where it
/// gets made, and deliberately:
///
///   * **Variation has to come from somewhere.** Evolution is variation plus
///     differential persistence. Phase 5h supplied the second. Without the
///     first, every daughter is its parent forever and nothing can happen.
///   * **Determinism survives untouched.** The engine has no global RNG.
///     Randomness is *derived* — `f(seed, domain, where, when)` — so a resumed
///     world draws the identical numbers, and the same seed gives the same
///     universe. What arrives here is stochasticity, not nondeterminism. The two
///     are routinely confused and are not the same thing.
///
/// # Fidelity is a consequence, not a setting
///
/// Each species is split by a binomial draw, so a species present in `n` copies
/// is transmitted with a relative error of `1/(2√n)`. Abundant molecules are
/// inherited almost exactly; rare ones are inherited badly or lost. **Copy
/// number is heredity's fidelity**, and it falls out of counting rather than
/// being declared — which is the same reason a real cell carries one genome and
/// not one molecule of each protein.
///
/// `p_a` is the probability a molecule lands in daughter A, which for a
/// symmetric pinch is the volume fraction and for the equal case is one half.
/// Nothing is created or destroyed: `a[i] + b[i] == counts[i]`, exactly, for
/// every species.
pub fn partition(counts: &[i128], p_a: f64, rng: &mut Rng) -> (Vec<i128>, Vec<i128>) {
    let p = p_a.clamp(0.0, 1.0);
    let mut a = Vec::with_capacity(counts.len());
    let mut b = Vec::with_capacity(counts.len());
    for &n in counts {
        let to_a = binomial(n.max(0), p, rng);
        a.push(to_a);
        b.push(n.max(0) - to_a);
    }
    (a, b)
}

/// How many of `n` land in A, each independently with probability `p`.
///
/// Two regimes, and the switch is stated rather than hidden. Up to
/// `EXACT_BELOW` the draw is exact: for `p = 1/2` that is a population count of
/// random bits, which is what a binomial *is*; otherwise it is a Bernoulli per
/// item. Above it the normal approximation is used, where the relative error it
/// commits is already below `1/√n < 1.6%` and falling — far under anything the
/// simulation can observe, since the quantity being approximated is itself the
/// noise.
fn binomial(n: i128, p: f64, rng: &mut Rng) -> i128 {
    const EXACT_BELOW: i128 = 1024;
    if n <= 0 {
        return 0;
    }
    if p <= 0.0 {
        return 0;
    }
    if p >= 1.0 {
        return n;
    }
    if n <= EXACT_BELOW {
        if (p - 0.5).abs() < 1e-12 {
            // A fair binomial is the number of ones in n random bits.
            let mut left = n;
            let mut k = 0i128;
            while left > 0 {
                let take = left.min(64);
                let word = rng.next_u64();
                let mask =
                    if take == 64 { u64::MAX } else { (1u64 << take) - 1 };
                k += (word & mask).count_ones() as i128;
                left -= take;
            }
            return k;
        }
        let mut k = 0i128;
        for _ in 0..n {
            if rng.chance(p) {
                k += 1;
            }
        }
        return k;
    }
    let mean = n as f64 * p;
    let sd = (n as f64 * p * (1.0 - p)).sqrt();
    let z = standard_normal(rng);
    let k = (mean + sd * z).round();
    (k as i128).clamp(0, n)
}

/// Box–Muller. Two uniforms in, one standard normal out.
fn standard_normal(rng: &mut Rng) -> f64 {
    // next_f64 may return exactly 0, and ln(0) is not a number we want.
    let u1 = rng.next_f64().max(f64::MIN_POSITIVE);
    let u2 = rng.next_f64();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// How much two compositions differ, as a number in `[0, 1]`.
///
/// The cosine distance between the count vectors: 0 for compositions with the
/// same proportions whatever their absolute size, 1 for compositions sharing no
/// species. Size-independent on purpose — a daughter is half its parent by
/// count and should not be scored as different for that reason alone. What is
/// inherited here is a *recipe*, not an amount.
pub fn composition_distance(x: &[i128], y: &[i128]) -> f64 {
    let mut dot = 0.0f64;
    let mut nx = 0.0f64;
    let mut ny = 0.0f64;
    for i in 0..x.len().max(y.len()) {
        let a = x.get(i).copied().unwrap_or(0) as f64;
        let b = y.get(i).copied().unwrap_or(0) as f64;
        dot += a * b;
        nx += a * a;
        ny += b * b;
    }
    if nx <= 0.0 || ny <= 0.0 {
        return 1.0;
    }
    (1.0 - dot / (nx.sqrt() * ny.sqrt())).clamp(0.0, 1.0)
}
