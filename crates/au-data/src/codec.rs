//! Binary serialization, hand-rolled.
//!
//! # Why not serde?
//!
//! Because the save format is not an implementation detail here — it is a
//! contract that has to survive a project measured in years. Between now and
//! Phase 13 every data structure in this repo will change shape, probably more
//! than once. If old saves cannot be migrated, then every schema change costs
//! us the world we were testing against, and we will start avoiding schema
//! changes to protect our saves. That is how a codebase ossifies.
//!
//! So: explicit little-endian encoding, an explicit schema version, and an
//! explicit place for migrations to live. It is a hundred lines and it means we
//! are never afraid of our own data.
//!
//! Little-endian everywhere so a save from an ARM laptop opens on an x86 server.

use std::fmt;

/// Bumped whenever the layout of anything written here changes.
/// Migrations key off this.
pub const FORMAT_VERSION: u16 = 3;

pub const MAGIC: &[u8; 4] = b"AUNI";

#[derive(Debug, PartialEq, Eq)]
pub enum CodecError {
    UnexpectedEof { needed: usize, had: usize },
    BadMagic,
    /// A file from a future version, or one we no longer know how to migrate.
    UnsupportedVersion(u16),
    ChecksumMismatch { expected: u64, found: u64 },
    /// A column id that no longer exists in the registry.
    UnknownColumn(u16),
    Invalid(&'static str),
}

impl fmt::Display for CodecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CodecError::UnexpectedEof { needed, had } => {
                write!(f, "truncated data: needed {} bytes, had {}", needed, had)
            }
            CodecError::BadMagic => write!(f, "not an Artificial Universe file"),
            CodecError::UnsupportedVersion(v) => write!(
                f,
                "save format v{} cannot be read by this build (current v{}); \
                 a migration is needed",
                v, FORMAT_VERSION
            ),
            CodecError::ChecksumMismatch { expected, found } => write!(
                f,
                "corrupt save: checksum {:#018x}, expected {:#018x}",
                found, expected
            ),
            CodecError::UnknownColumn(id) => {
                write!(f, "column {:#06x} is not in the registry", id)
            }
            CodecError::Invalid(m) => write!(f, "invalid data: {}", m),
        }
    }
}

impl std::error::Error for CodecError {}

#[derive(Default)]
pub struct Writer {
    buf: Vec<u8>,
}

impl Writer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn u8(&mut self, v: u8) -> &mut Self {
        self.buf.push(v);
        self
    }
    pub fn u16(&mut self, v: u16) -> &mut Self {
        self.buf.extend_from_slice(&v.to_le_bytes());
        self
    }
    pub fn u32(&mut self, v: u32) -> &mut Self {
        self.buf.extend_from_slice(&v.to_le_bytes());
        self
    }
    pub fn i32(&mut self, v: i32) -> &mut Self {
        self.buf.extend_from_slice(&v.to_le_bytes());
        self
    }
    pub fn u64(&mut self, v: u64) -> &mut Self {
        self.buf.extend_from_slice(&v.to_le_bytes());
        self
    }
    /// 16 bytes. The exact conserved quantities (energy, mass) live at this
    /// width — see `au_physics::quantity` for why they are not floats.
    pub fn i128(&mut self, v: i128) -> &mut Self {
        self.buf.extend_from_slice(&v.to_le_bytes());
        self
    }
    pub fn f32(&mut self, v: f32) -> &mut Self {
        self.buf.extend_from_slice(&v.to_bits().to_le_bytes());
        self
    }
    /// Floats are written by their bit pattern, never by a decimal
    /// representation. A round-trip must be bit-exact or determinism is gone.
    pub fn f64(&mut self, v: f64) -> &mut Self {
        self.buf.extend_from_slice(&v.to_bits().to_le_bytes());
        self
    }
    pub fn bytes(&mut self, b: &[u8]) -> &mut Self {
        self.buf.extend_from_slice(b);
        self
    }

    pub fn len(&self) -> usize {
        self.buf.len()
    }
    pub fn is_empty(&self) -> bool {
        self.buf.is_empty()
    }
    pub fn as_slice(&self) -> &[u8] {
        &self.buf
    }
    pub fn finish(self) -> Vec<u8> {
        self.buf
    }
}

pub struct Reader<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    pub fn new(buf: &'a [u8]) -> Self {
        Reader { buf, pos: 0 }
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], CodecError> {
        if self.pos + n > self.buf.len() {
            return Err(CodecError::UnexpectedEof {
                needed: n,
                had: self.buf.len() - self.pos,
            });
        }
        let s = &self.buf[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    pub fn u8(&mut self) -> Result<u8, CodecError> {
        Ok(self.take(1)?[0])
    }
    pub fn u16(&mut self) -> Result<u16, CodecError> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }
    pub fn u32(&mut self) -> Result<u32, CodecError> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    pub fn i32(&mut self) -> Result<i32, CodecError> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    pub fn u64(&mut self) -> Result<u64, CodecError> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    pub fn i128(&mut self) -> Result<i128, CodecError> {
        Ok(i128::from_le_bytes(self.take(16)?.try_into().unwrap()))
    }
    pub fn f32(&mut self) -> Result<f32, CodecError> {
        Ok(f32::from_bits(u32::from_le_bytes(self.take(4)?.try_into().unwrap())))
    }
    pub fn f64(&mut self) -> Result<f64, CodecError> {
        Ok(f64::from_bits(u64::from_le_bytes(self.take(8)?.try_into().unwrap())))
    }
    pub fn bytes(&mut self, n: usize) -> Result<&'a [u8], CodecError> {
        self.take(n)
    }

    pub fn remaining(&self) -> usize {
        self.buf.len() - self.pos
    }
    pub fn pos(&self) -> usize {
        self.pos
    }
}
