//! Deep time, exactly.
//!
//! # Why not `f64` seconds?
//!
//! This simulation must run for billions of simulated years while still being
//! able to resolve a chemical reaction. An `f64` has 53 bits of mantissa; at
//! 4.5e9 years (~1.4e17 s) it can no longer represent a single second, let
//! alone a fraction of one. Worse, accumulating `now += dt` in floating point
//! drifts, and drift means two runs of the same seed diverge. Determinism is
//! the property everything else in this project rests on.
//!
//! So time is an **exact 128-bit fixed-point value**: `u64` seconds plus a
//! `u64` fraction in units of 1/2^64 of a second. Every advance is integer
//! addition with carry. Bit-exact, replayable, and good for ~5.8e11 years with
//! sub-attosecond resolution. It cannot drift, because there is nothing to
//! round.
//!
//! # Two different "speeds" — do not confuse them
//!
//! * [`SimClock::scale`] — **simulated seconds per tick.** This is simulation
//!   state. Changing it changes what the universe computes.
//! * Wall-clock speed (ticks executed per real-world second) — **presentation
//!   only.** It lives in the runner/renderer and is invisible to the sim.
//!
//! The player pressing fast-forward must never alter history. That is why the
//! player's speed control is not in this file, and must never be.

/// Seconds in a Julian year (365.25 d). Fixed by definition; not a tunable.
pub const SECONDS_PER_YEAR: f64 = 31_557_600.0;

/// 2^64, as an f64. Used only when converting *from* human-authored config.
const FRAC_SCALE: f64 = 18_446_744_073_709_551_616.0;

/// A monotonically increasing simulation step counter.
///
/// The tick, not the wall clock and not `SimInstant`, is the canonical index
/// of "when" something happened. Scheduling is expressed in ticks.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct Tick(pub u64);

impl Tick {
    #[inline]
    pub fn next(self) -> Tick {
        Tick(self.0 + 1)
    }
}

/// An exact span of simulated time. 128-bit fixed point.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct SimDuration {
    pub secs: u64,
    /// Fraction of a second, in units of 1/2^64 s.
    pub frac: u64,
}

impl SimDuration {
    pub const ZERO: SimDuration = SimDuration { secs: 0, frac: 0 };

    pub const fn from_secs(secs: u64) -> Self {
        SimDuration { secs, frac: 0 }
    }

    /// Convert from a human-authored `f64` (config files, tools).
    ///
    /// This is the *only* lossy door into the time system, and it is one-way:
    /// once a duration exists it is exact forever. Never call this inside the
    /// tick loop.
    pub fn from_secs_f64(secs: f64) -> Self {
        assert!(secs.is_finite() && secs >= 0.0, "SimDuration must be finite and non-negative");
        let whole = secs.trunc();
        assert!(whole < u64::MAX as f64, "SimDuration overflow");
        SimDuration {
            secs: whole as u64,
            frac: (secs.fract() * FRAC_SCALE) as u64,
        }
    }

    pub fn from_years(years: f64) -> Self {
        Self::from_secs_f64(years * SECONDS_PER_YEAR)
    }

    /// Lossy. For display and for coarse physics that has already accepted
    /// float error. Never feed the result back into the clock.
    pub fn as_secs_f64(self) -> f64 {
        self.secs as f64 + (self.frac as f64) / FRAC_SCALE
    }

    pub fn as_years_f64(self) -> f64 {
        self.as_secs_f64() / SECONDS_PER_YEAR
    }

    /// Exact addition with carry. Panics on overflow rather than wrapping:
    /// a silently wrapped clock would corrupt every downstream layer.
    pub fn add(self, other: SimDuration) -> SimDuration {
        let (frac, carry) = self.frac.overflowing_add(other.frac);
        let secs = self
            .secs
            .checked_add(other.secs)
            .and_then(|s| s.checked_add(carry as u64))
            .expect("simulation time overflow (>5.8e11 years)");
        SimDuration { secs, frac }
    }

    /// Exact multiplication by an integer count of ticks.
    ///
    /// Used to fast-forward the clock across a skipped span without looping.
    pub fn mul_u64(self, n: u64) -> SimDuration {
        let total_frac = (self.frac as u128) * (n as u128);
        let carry = (total_frac >> 64) as u64;
        let frac = total_frac as u64;
        let secs = self
            .secs
            .checked_mul(n)
            .and_then(|s| s.checked_add(carry))
            .expect("simulation time overflow");
        SimDuration { secs, frac }
    }
}

/// An exact point in simulated time, measured from the birth of the universe.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct SimInstant {
    pub secs: u64,
    pub frac: u64,
}

impl SimInstant {
    pub const ORIGIN: SimInstant = SimInstant { secs: 0, frac: 0 };

    pub fn advanced_by(self, d: SimDuration) -> SimInstant {
        let sum = SimDuration { secs: self.secs, frac: self.frac }.add(d);
        SimInstant { secs: sum.secs, frac: sum.frac }
    }

    pub fn as_secs_f64(self) -> f64 {
        self.secs as f64 + (self.frac as f64) / FRAC_SCALE
    }

    pub fn as_years_f64(self) -> f64 {
        self.as_secs_f64() / SECONDS_PER_YEAR
    }
}

/// The universal simulation clock (ARCHITECTURE §3).
///
/// Holds the tick counter, the exact instant, and the current time *scale*
/// (simulated seconds per tick). The scale is simulation state: it is
/// snapshotted, hashed, and replayed. It exists so that the engine can run
/// chemistry at seconds-per-tick and deep evolution at millennia-per-tick
/// without changing any system's code — only its scheduling period and the
/// clock's scale.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct SimClock {
    tick: Tick,
    now: SimInstant,
    scale: SimDuration,
}

impl SimClock {
    pub fn new(scale: SimDuration) -> Self {
        SimClock { tick: Tick(0), now: SimInstant::ORIGIN, scale }
    }

    #[inline]
    pub fn tick(&self) -> Tick {
        self.tick
    }

    #[inline]
    pub fn now(&self) -> SimInstant {
        self.now
    }

    #[inline]
    pub fn scale(&self) -> SimDuration {
        self.scale
    }

    /// Advance exactly one tick.
    #[inline]
    pub fn advance(&mut self) {
        self.tick = self.tick.next();
        self.now = self.now.advanced_by(self.scale);
    }

    /// Change simulated seconds-per-tick.
    ///
    /// This is a **simulation event**, not a UI action. It alters the outcome
    /// of the universe and is recorded as such. Player fast-forward does not
    /// call this.
    pub fn set_scale(&mut self, scale: SimDuration) {
        self.scale = scale;
    }

    /// Restore a clock from a snapshot without replaying its history.
    pub fn restore(tick: Tick, now: SimInstant, scale: SimDuration) -> Self {
        SimClock { tick, now, scale }
    }
}
