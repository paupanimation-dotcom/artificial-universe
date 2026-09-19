//! # au-data — the data architecture
//!
//! Answers ARCHITECTURE §19 (data architecture) and §20 (performance) with
//! three ideas:
//!
//! 1. **Structure of arrays.** Layers stream one field at a time through cache.
//! 2. **Dirty tracking + procedural reconstruction.** Untouched world is not
//!    stored; it is a pure function of the seed, and can be dropped and rebuilt
//!    bit-identically. Only *changes* cost memory.
//! 3. **Level of detail.** Chunks are simulated fully, as populations, or as
//!    statistics. The player zooms from one organism to a planet; the engine
//!    pays only where it must.
//!
//! Together these decide whether "the universe must be enormous" is a design
//! goal or a fantasy.

pub mod chunk;
pub mod codec;
pub mod column;
pub mod snapshot;

pub use chunk::{Chunk, ChunkCoord, ChunkGenerator, ChunkStore, EmptyGenerator, Lod};
pub use codec::{CodecError, Reader, Writer, FORMAT_VERSION};
pub use column::{Column, ColumnId, ColumnRegistry, ColumnSet, Field, VecColumn};
