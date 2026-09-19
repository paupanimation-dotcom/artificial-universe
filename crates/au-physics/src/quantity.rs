//! Conserved quantities, as exact integers.
//!
//! # Why integers, and why this is the most important decision in Phase 2
//!
//! **Evolution is an adversarial optimizer, and its search space includes the
//! bugs in your physics engine.**
//!
//! If energy is not conserved — if there is *any* path by which a cell can end
//! up with more energy than it was given — then sooner or later something will
//! evolve to walk that path. Not because it is clever, but because natural
//! selection is an exhaustive search and a free-energy exploit is the single
//! most rewarding thing in any possible fitness landscape. This is not a
//! hypothetical failure mode; it is the classic way artificial-life projects
//! die. You end up with organisms whose entire strategy is a floating-point
//! rounding error, and no amount of ecology on top will fix it.
//!
//! Floating point cannot give us conservation. `a -= f; b += f;` in `f64` does
//! not preserve `a + b`: both additions round. Over 10^9 ticks that leak is not
//! small, and it is not random — it is *systematically exploitable*.
//!
//! So the conserved quantities are integers, and transport is **symmetric
//! flux**: the same integer is subtracted from one cell and added to the other.
//! Conservation then holds *to the bit*, forever, as an identity rather than an
//! approximation. There is a test that runs a hundred thousand ticks and asserts
//! the total is exactly, precisely unchanged.
//!
//! # Where floats are still allowed
//!
//! Temperature, density, conductivity, pressure — all **derived**, all
//! recomputed from the exact integers every tick, none ever stored. Float error
//! cannot accumulate if you never accumulate floats. The error in a derived
//! value is bounded by a single operation and dies at the end of the tick.
//!
//! # Units and range
//!
//! * Energy: microjoules. `i128` ⇒ ±1.7 × 10^32 J. (Earth's entire thermal
//!   energy content is ~10^30 J, so a planet fits with room to spare, while a
//!   micromole of chemistry still resolves.)
//! * Mass: micrograms. `i128` ⇒ ±1.7 × 10^29 kg. (Earth is 6 × 10^24 kg.)
//!
//! A single quantum is far below anything the simulation will ever need to
//! resolve, and far above the point where `f64` would have started lying.

use au_core::hash::{Hasher, WorldHash};

/// Microjoules per joule.
/// Energy quantum: **picojoules** per joule (10¹²).
///
/// Raised from microjoules for the same reason the mass quantum was raised from
/// micrograms, and discovered the same way — by watching a cell fail to exist.
/// With mass fixed, a cubic micron of rock at 1500 K still held *zero* energy:
/// `m·c·ΔT` is 3.75×10⁻⁹ J there, four thousandths of one microjoule. No energy
/// means no derived temperature, no phase, an immobile cell, and a membrane that
/// cannot form. The foundation had two floors, not one.
///
/// # Why not finer, which is the interesting part
///
/// The obvious move is to go as fine as `i128` allows. It is wrong, and the
/// reason is a second ceiling nobody declares: `as_joules` divides through
/// `f64`, and an `f64` only holds integers exactly up to 2⁵³ ≈ 9×10¹⁵. Past
/// that, the low bits of a stored energy are gone.
///
/// That matters because the quantity we usually care about is a *difference*
/// against a large baseline. A convection cell sits at 300 K and is driven by a
/// temperature perturbation that may be a millionth of a kelvin; the absolute
/// energy is enormous compared with the difference that does the work. Lose the
/// low bits and every cell reads the same temperature, buoyancy becomes exactly
/// zero, and the fluid simply stops convecting. Femtojoules were tried first and
/// did precisely that — the momentum books still balanced perfectly, because
/// nothing was leaking; there was just nothing to push with.
///
/// So the quantum is pinned between two hard limits:
///
/// ```text
///     quantum        fluid-tank cell    micron cell
///     microjoule     1.2e9   ok         0.004     rounds to nothing
///     nanojoule      1.2e12  ok         3.75      too coarse to be a phase
///     picojoule      1.2e15  ok         3.75e3    ok
///     femtojoule     1.2e18  PAST 2^53  3.75e6    ok, but the fluid dies
/// ```
///
/// Picojoules is the only value that clears both ends, and it clears the upper
/// one by about sevenfold. The old microjoule value happened to suit the tests
/// that existed; nothing chose it for that reason.
///
/// The deeper fix, when this becomes binding again, is to store energy
/// *relative to a reference* instead of absolutely, which removes the
/// cancellation entirely. That changes the conservation invariants and deserves
/// its own phase.
pub const PJ_PER_J: i128 = 1_000_000_000_000;
/// Micrograms per kilogram.
/// Mass quantum: **attograms** per kilogram (10²¹).
///
/// # Why this is not micrograms any more
///
/// It was, and that turned out to cap the engine's spatial resolution at about
/// a tenth of a millimetre — a limit discovered from the far end, when
/// membranes refused to work. A cubic micron of rock weighs 2.5×10⁻¹⁵ kg, which
/// is a *millionth* of one microgram. It rounded to zero, so the host had no
/// mass, so its phase could not be derived, so the mobility mask called the cell
/// immobile and nothing moved. The chemistry was fine; the foundation could not
/// hold a cell.
///
/// ```text
///     cell        side     micrograms   attograms
///     10⁻¹⁸ m³    1 µm     0.0000025    2.5 × 10⁶
///     10⁻¹⁵ m³    10 µm    0.0025       2.5 × 10⁹
///     10⁻⁹  m³    1 mm     2500         2.5 × 10¹⁵
/// ```
///
/// # Choosing the quantum
///
/// An `i128` spans about 1.7×10³⁸, and that budget has to cover both ends of the
/// range the engine actually uses:
///
///   * the smallest thing worth resolving — a bacterium-scale cell of water at
///     10⁻¹⁵ kg, which should be many quanta rather than one;
///   * the largest single cell — planetary terrain at ~3×10⁶ kg for a 10 m cell,
///     or ~3×10¹² kg if someone builds a kilometre grid.
///
/// Attograms give a bacterium 10⁶ quanta and a ceiling of 1.7×10¹⁷ kg per cell —
/// room for grid cells up to roughly 30 km on a side. Finer (yoctograms) would
/// resolve individual molecules but overflow on planetary cells; coarser
/// (picograms) puts us back where we started. Note that molecular masses are
/// *not* resolved here and do not need to be: chemistry counts molecules as
/// integers of its own, and this quantity is the bulk mass of the host.
pub const AG_PER_KG: i128 = 1_000_000_000_000_000_000_000;

/// Exact energy. Microjoules.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct Energy(pub i128);

impl Energy {
    pub const ZERO: Energy = Energy(0);

    /// The one lossy door in. Rounds to the nearest microjoule, deterministically
    /// (`f64::round` is round-half-away-from-zero on every IEEE-754 host).
    /// Saturates rather than wrapping — a wrapped energy would be a free-energy
    /// exploit handed to evolution on a plate.
    #[inline]
    pub fn from_joules(j: f64) -> Energy {
        debug_assert!(j.is_finite(), "non-finite energy: {}", j);
        Energy((j * PJ_PER_J as f64).round() as i128)
    }

    #[inline]
    pub fn as_joules(self) -> f64 {
        self.0 as f64 / PJ_PER_J as f64
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.0 == 0
    }
}

impl std::ops::Add for Energy {
    type Output = Energy;
    #[inline]
    fn add(self, o: Energy) -> Energy {
        Energy(self.0.checked_add(o.0).expect("energy overflow"))
    }
}

impl std::ops::Sub for Energy {
    type Output = Energy;
    #[inline]
    fn sub(self, o: Energy) -> Energy {
        Energy(self.0.checked_sub(o.0).expect("energy underflow"))
    }
}

impl std::ops::AddAssign for Energy {
    #[inline]
    fn add_assign(&mut self, o: Energy) {
        *self = *self + o;
    }
}

impl std::ops::SubAssign for Energy {
    #[inline]
    fn sub_assign(&mut self, o: Energy) {
        *self = *self - o;
    }
}

impl WorldHash for Energy {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_bytes(&self.0.to_le_bytes());
    }
}

/// Momentum, in units of 10⁻¹⁵ kg·m·s⁻¹ (picogram-metres per second).
///
/// Exact, for the same reason energy is exact — and the reason is worth stating
/// separately, because it is not the same reason.
///
/// Energy must be conserved or something evolves to eat the leak. **Momentum
/// must be conserved or something evolves to *swim* on the leak.** A fluid that
/// can gain momentum from nothing is a reactionless drive, and a reactionless
/// drive is free propulsion — which, to a selection process, is indistinguishable
/// from free food.
///
/// Symmetric flux again: whatever momentum one parcel takes from another, the
/// other loses exactly. Newton's third law, as an integer identity rather than an
/// aspiration.
///
/// Range: ±1.7 × 10²³ kg·m·s⁻¹ per cell. Resolution 10⁻¹⁵ — fine enough to
/// resolve mantle creep at 10⁻¹⁰ m·s⁻¹.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct Momentum(pub i128);

/// Momentum units per kg·m·s⁻¹.
pub const P_PER_KG_M_S: i128 = 1_000_000_000_000_000;

impl Momentum {
    pub const ZERO: Momentum = Momentum(0);

    #[inline]
    pub fn from_si(p: f64) -> Momentum {
        debug_assert!(p.is_finite(), "non-finite momentum: {}", p);
        Momentum((p * P_PER_KG_M_S as f64).round() as i128)
    }

    #[inline]
    pub fn as_si(self) -> f64 {
        self.0 as f64 / P_PER_KG_M_S as f64
    }
}

impl std::ops::AddAssign for Momentum {
    #[inline]
    fn add_assign(&mut self, o: Momentum) {
        self.0 = self.0.checked_add(o.0).expect("momentum overflow");
    }
}

impl WorldHash for Momentum {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_bytes(&self.0.to_le_bytes());
    }
}

/// Exact mass. Micrograms.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct Mass(pub i128);

impl Mass {
    pub const ZERO: Mass = Mass(0);

    #[inline]
    pub fn from_kg(kg: f64) -> Mass {
        debug_assert!(kg.is_finite() && kg >= 0.0, "bad mass: {}", kg);
        Mass((kg * AG_PER_KG as f64).round() as i128)
    }

    #[inline]
    pub fn as_kg(self) -> f64 {
        self.0 as f64 / AG_PER_KG as f64
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.0 == 0
    }
}

impl WorldHash for Mass {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_bytes(&self.0.to_le_bytes());
    }
}

/// The books.
///
/// A closed system conserves energy. An *open* one — a planet lit by a star and
/// radiating into the night — does not, and must not: the whole reason life is
/// possible is that energy flows *through*. So we cannot simply assert that the
/// total never changes. We assert something stronger and more useful:
///
/// ```text
/// E_total(now) - E_total(start)  ==  energy_in - energy_out     EXACTLY
/// ```
///
/// Every joule that crosses a boundary is written down. Nothing appears, nothing
/// vanishes, and if either ever happens we know on the very next tick instead of
/// discovering it in Phase 9 when something evolves to live on it.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct Ledger {
    pub energy_in: Energy,
    pub energy_out: Energy,
    pub mass_in: Mass,
    pub mass_out: Mass,

    /// Total external impulse delivered to the fluid, per axis.
    ///
    /// Advection and viscosity are symmetric exchanges and contribute exactly
    /// zero here — they only move momentum around. What lands in this ledger is
    /// everything that came from *outside*: buoyancy, the pressure force on the
    /// walls, friction against a no-slip surface, momentum absorbed by a solid.
    ///
    /// Which gives the identity that the fluid solver lives or dies by:
    ///
    /// ```text
    /// Σ momentum(now) − Σ momentum(start)  ==  impulse     EXACTLY
    /// ```
    ///
    /// Break it and the fluid can push on nothing. There is a test.
    pub impulse: [Momentum; 3],
}

impl Ledger {
    /// Net energy the world should have gained since the beginning.
    pub fn net_energy(&self) -> Energy {
        self.energy_in - self.energy_out
    }
}

impl WorldHash for Ledger {
    fn hash_into(&self, h: &mut Hasher) {
        self.energy_in.hash_into(h);
        self.energy_out.hash_into(h);
        self.mass_in.hash_into(h);
        self.mass_out.hash_into(h);
        for i in self.impulse {
            i.hash_into(h);
        }
    }
}
