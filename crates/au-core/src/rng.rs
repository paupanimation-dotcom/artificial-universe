//! Deterministic, splittable randomness.
//!
//! # The rule
//!
//! **There is no global RNG. Ever.**
//!
//! A single shared random stream is the fastest way to destroy this project.
//! The moment two systems draw from the same generator, the *order* in which
//! they happen to run becomes part of the result — and the moment we run
//! chunks in parallel (which we will, by Phase 4), that order stops being
//! deterministic. The universe would stop being replayable, family trees would
//! stop reconstructing, and the bug would look like "evolution is slightly
//! different sometimes", which is the worst possible bug to have to find.
//!
//! Instead, randomness is **derived**, not drawn:
//!
//! ```text
//! Rng = f(world_seed, domain, key_a, key_b)
//! ```
//!
//! Every system, at every chunk, on every tick, derives its own private stream
//! from pure inputs. Nothing is shared, so nothing needs locking, so ordering
//! is irrelevant, so parallelism is free and determinism is preserved.
//!
//! A mutation in the year 4,000,000 is reproducible from the seed alone.

/// SplitMix64 — the seed expander. Turns a small integer into well-mixed bits.
#[inline]
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Namespace tag for a random stream.
///
/// Two different domains asking for "the randomness at chunk C on tick T" must
/// get independent streams. Add a constant here whenever a new system needs
/// randomness; never reuse another system's domain.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub struct Domain(pub u64);

impl Domain {
    pub const WORLD_GEN: Domain = Domain(0x0001);
    pub const CHUNK_GEN: Domain = Domain(0x0002);
    // Reserved for future layers. Kept here so collisions are visible at a glance.
    pub const PHYSICS: Domain = Domain(0x0100);
    pub const CHEMISTRY: Domain = Domain(0x0200);
    pub const PLANET: Domain = Domain(0x0400);
    pub const MUTATION: Domain = Domain(0x0300);
    pub const RECOMBINATION: Domain = Domain(0x0301);
    pub const BEHAVIOR: Domain = Domain(0x0400);
    /// For tests only.
    pub const TEST: Domain = Domain(0xFFFF);
}

/// xoshiro256++ — fast, well-distributed, 256 bits of state.
///
/// Not cryptographic. Does not need to be. It needs to be *identical on every
/// machine*, which is why it is written out here in integer ops rather than
/// pulled from a crate whose next minor version might change an algorithm.
#[derive(Clone, Debug)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    /// Derive a private stream from pure inputs.
    ///
    /// This is the only way to make an `Rng`. Same inputs, same stream, on any
    /// machine, forever.
    pub fn derive(world_seed: u64, domain: Domain, key_a: u64, key_b: u64) -> Rng {
        // Mix all four inputs before expansion so that adjacent keys (e.g.
        // neighbouring chunks, consecutive ticks) produce completely unrelated
        // streams rather than correlated ones.
        let mut z = world_seed
            ^ domain.0.wrapping_mul(0xD6E8_FEB8_6659_FD93)
            ^ key_a.rotate_left(17).wrapping_mul(0xA076_1D64_78BD_642F)
            ^ key_b.rotate_left(41).wrapping_mul(0xE703_7ED1_A0B4_28DB);
        let s = [
            splitmix64(&mut z),
            splitmix64(&mut z),
            splitmix64(&mut z),
            splitmix64(&mut z),
        ];
        Rng { s }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    #[inline]
    pub fn next_u32(&mut self) -> u32 {
        (self.next_u64() >> 32) as u32
    }

    /// Uniform in [0, 1). Exactly 53 bits of mantissa — the most an f64 holds.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        // Take the top 53 bits; scale by 2^-53. Deterministic on any IEEE-754 host.
        ((self.next_u64() >> 11) as f64) * (1.0 / 9_007_199_254_740_992.0)
    }

    /// Uniform in [0, n). Rejection-sampled, so it is unbiased *and* its bias
    /// does not depend on the platform.
    pub fn below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "Rng::below requires n > 0");
        let zone = u64::MAX - (u64::MAX % n);
        loop {
            let v = self.next_u64();
            if v < zone {
                return v % n;
            }
        }
    }

    /// True with probability `p`.
    #[inline]
    pub fn chance(&mut self, p: f64) -> bool {
        self.next_f64() < p
    }
}
