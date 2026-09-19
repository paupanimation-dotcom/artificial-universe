//! World hashing — the instrument that makes determinism testable.
//!
//! An emergent simulation is, by construction, a system whose behaviour nobody
//! can predict by reading the code. That makes normal debugging weak: you
//! cannot assert that a wolf appears, because there are no wolves and there is
//! no "should".
//!
//! What you *can* assert is that the same seed produces the same universe. So
//! every piece of simulation state folds into a single 64-bit digest. Two runs
//! agree, or they do not. When they diverge, we bisect by tick and layer until
//! we find the exact system that stopped being reproducible.
//!
//! This is the single most valuable diagnostic in the whole project, and it is
//! ~40 lines. Written on day one, before there is anything to hash, because by
//! Phase 8 the causal chain will be 10^9 ticks long and retrofitting it will be
//! impossible.
//!
//! FNV-1a: order-sensitive by design. If chunk iteration order changes, the
//! hash changes — and we *want* to be told about that, because a world whose
//! iteration order is unstable is a world that cannot be replayed.

#[derive(Clone, Copy, Debug)]
pub struct Hasher(u64);

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

impl Default for Hasher {
    fn default() -> Self {
        Hasher(FNV_OFFSET)
    }
}

impl Hasher {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn write_u8(&mut self, v: u8) {
        self.0 ^= v as u64;
        self.0 = self.0.wrapping_mul(FNV_PRIME);
    }

    #[inline]
    pub fn write_bytes(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.write_u8(b);
        }
    }

    #[inline]
    pub fn write_u16(&mut self, v: u16) {
        self.write_bytes(&v.to_le_bytes());
    }

    #[inline]
    pub fn write_u32(&mut self, v: u32) {
        self.write_bytes(&v.to_le_bytes());
    }

    #[inline]
    pub fn write_u64(&mut self, v: u64) {
        self.write_bytes(&v.to_le_bytes());
    }

    #[inline]
    pub fn write_i32(&mut self, v: i32) {
        self.write_bytes(&v.to_le_bytes());
    }

    /// Hash a float by its **bits**, not its value.
    ///
    /// Deliberate: `0.1 + 0.2` and `0.3` must hash differently, because they
    /// *are* different, and a physics change that flips the last mantissa bit
    /// is exactly the kind of drift we are hunting. Normalise NaN so that the
    /// many NaN bit-patterns don't create phantom divergence.
    #[inline]
    pub fn write_f64(&mut self, v: f64) {
        let bits = if v.is_nan() { f64::NAN.to_bits() } else { v.to_bits() };
        self.write_u64(bits);
    }

    #[inline]
    pub fn finish(self) -> u64 {
        self.0
    }
}

/// Anything that contributes to the identity of the world.
pub trait WorldHash {
    fn hash_into(&self, h: &mut Hasher);

    fn world_hash(&self) -> u64
    where
        Self: Sized,
    {
        let mut h = Hasher::new();
        self.hash_into(&mut h);
        h.finish()
    }
}
