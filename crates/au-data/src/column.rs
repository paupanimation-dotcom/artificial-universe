//! Structure-of-arrays storage.
//!
//! # Why SoA and not "a Vec of Organism structs"
//!
//! Because of the numbers. This project wants millions of entities and a tick
//! rate that makes deep time reachable. A physics pass that touches only
//! `temperature` should stream *only* temperature through the cache — not drag
//! along genome pointers, brain state and family-tree links on every cache line.
//! The difference between array-of-structs and structure-of-arrays here is not
//! a micro-optimisation; it is one to two orders of magnitude, and it decides
//! whether Phase 9 is a simulation or a slideshow.
//!
//! It is also what makes SIMD and multi-threading possible later without a
//! rewrite. Getting this wrong is the kind of mistake that is invisible for two
//! years and then unfixable.
//!
//! # The erasure
//!
//! Chunks must hold *whatever* columns the layers above define — and `au-data`
//! must not know what those are. So a column is a `dyn Column` trait object
//! with a numeric [`ColumnId`], and decoding goes through a [`ColumnRegistry`]
//! that maps ids back to constructors. Typed access is recovered with `Any`
//! downcasting: safe, no `unsafe`, and the cost is one vtable hop per *column*,
//! not per element — the hot loop still iterates a flat `&[T]`.

use au_core::hash::{Hasher, WorldHash};
use std::any::Any;
use std::collections::BTreeMap;

use crate::codec::{CodecError, Reader, Writer};

/// Numeric identity of a column. Assigned by the layer that owns it.
///
/// Ids are permanent once shipped: a save file refers to columns by number, so
/// renumbering a column silently reinterprets old worlds. Reserve ranges per
/// layer (0x1xxx physics, 0x2xxx chemistry, ...) so two layers never collide.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub struct ColumnId(pub u16);

/// A field that can live in a column.
///
/// Bit-exact round-trip is mandatory: `decode(encode(x)) == x`, always.
pub trait Field: Copy + 'static {
    const ID: ColumnId;
    fn encode(&self, w: &mut Writer);
    fn decode(r: &mut Reader) -> Result<Self, CodecError>;
    fn hash_into(&self, h: &mut Hasher);
}

// Blanket implementations for the primitives every layer will want.
impl Field for f32 {
    const ID: ColumnId = ColumnId(0);
    fn encode(&self, w: &mut Writer) {
        w.f32(*self);
    }
    fn decode(r: &mut Reader) -> Result<Self, CodecError> {
        r.f32()
    }
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u32(self.to_bits());
    }
}

impl Field for f64 {
    const ID: ColumnId = ColumnId(0);
    fn encode(&self, w: &mut Writer) {
        w.f64(*self);
    }
    fn decode(r: &mut Reader) -> Result<Self, CodecError> {
        r.f64()
    }
    fn hash_into(&self, h: &mut Hasher) {
        h.write_f64(*self);
    }
}

impl Field for u32 {
    const ID: ColumnId = ColumnId(0);
    fn encode(&self, w: &mut Writer) {
        w.u32(*self);
    }
    fn decode(r: &mut Reader) -> Result<Self, CodecError> {
        r.u32()
    }
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u32(*self);
    }
}

/// The width the conserved quantities use. Sixteen bytes per cell per field is
/// not cheap — but the alternative is an `f64` that leaks energy, and a leak is
/// an exploit that something will eventually evolve to live on.
impl Field for i128 {
    const ID: ColumnId = ColumnId(0);
    fn encode(&self, w: &mut Writer) {
        w.i128(*self);
    }
    fn decode(r: &mut Reader) -> Result<Self, CodecError> {
        r.i128()
    }
    fn hash_into(&self, h: &mut Hasher) {
        h.write_bytes(&self.to_le_bytes());
    }
}

impl Field for u16 {
    const ID: ColumnId = ColumnId(0);
    fn encode(&self, w: &mut Writer) {
        w.u16(*self);
    }
    fn decode(r: &mut Reader) -> Result<Self, CodecError> {
        r.u16()
    }
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u16(*self);
    }
}

impl Field for u64 {
    const ID: ColumnId = ColumnId(0);
    fn encode(&self, w: &mut Writer) {
        w.u64(*self);
    }
    fn decode(r: &mut Reader) -> Result<Self, CodecError> {
        r.u64()
    }
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u64(*self);
    }
}

/// Type-erased column, so `au-data` can store fields it has never heard of.
pub trait Column: Any {
    fn id(&self) -> ColumnId;
    fn len(&self) -> usize;
    fn encode(&self, w: &mut Writer);
    fn hash_into(&self, h: &mut Hasher);
    fn as_any(&self) -> &dyn Any;
    fn as_any_mut(&mut self) -> &mut dyn Any;
    fn clone_box(&self) -> Box<dyn Column>;
}

/// A dense array of one field. The thing the hot loops actually touch.
#[derive(Clone, Debug)]
pub struct VecColumn<T: Field> {
    id: ColumnId,
    data: Vec<T>,
}

impl<T: Field> VecColumn<T> {
    /// `id` is explicit rather than taken from `T::ID` so that two columns of
    /// the same primitive type (say, `temperature: f64` and `pressure: f64`)
    /// remain distinct entities in storage and in save files.
    pub fn new(id: ColumnId) -> Self {
        VecColumn { id, data: Vec::new() }
    }

    pub fn with_data(id: ColumnId, data: Vec<T>) -> Self {
        VecColumn { id, data }
    }

    pub fn as_slice(&self) -> &[T] {
        &self.data
    }
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }
    pub fn push(&mut self, v: T) {
        self.data.push(v);
    }
}

impl<T: Field + std::fmt::Debug> Column for VecColumn<T> {
    fn id(&self) -> ColumnId {
        self.id
    }
    fn len(&self) -> usize {
        self.data.len()
    }
    fn encode(&self, w: &mut Writer) {
        w.u32(self.data.len() as u32);
        for v in &self.data {
            v.encode(w);
        }
    }
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u16(self.id.0);
        h.write_u32(self.data.len() as u32);
        for v in &self.data {
            v.hash_into(h);
        }
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
    fn clone_box(&self) -> Box<dyn Column> {
        Box::new(self.clone())
    }
}

/// Rebuilds a column from bytes. One per registered column id.
type DecodeFn = fn(ColumnId, &mut Reader) -> Result<Box<dyn Column>, CodecError>;

/// Maps column ids back to decoders.
///
/// Every layer registers its columns here at startup. A save file naming a
/// column nobody registered is an error, not a silent skip — because silently
/// dropping a column would load a world that *looks* fine and is subtly wrong,
/// which is the worst failure mode this project has.
#[derive(Default, Clone)]
pub struct ColumnRegistry {
    decoders: BTreeMap<ColumnId, DecodeFn>,
}

impl ColumnRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn register<T: Field + std::fmt::Debug>(&mut self, id: ColumnId) {
        fn decode<T: Field + std::fmt::Debug>(
            id: ColumnId,
            r: &mut Reader,
        ) -> Result<Box<dyn Column>, CodecError> {
            let n = r.u32()? as usize;
            let mut data = Vec::with_capacity(n);
            for _ in 0..n {
                data.push(T::decode(r)?);
            }
            Ok(Box::new(VecColumn::with_data(id, data)))
        }
        self.decoders.insert(id, decode::<T>);
    }

    pub fn decode(&self, id: ColumnId, r: &mut Reader) -> Result<Box<dyn Column>, CodecError> {
        let f = self.decoders.get(&id).ok_or(CodecError::UnknownColumn(id.0))?;
        f(id, r)
    }

    pub fn contains(&self, id: ColumnId) -> bool {
        self.decoders.contains_key(&id)
    }

    pub fn len(&self) -> usize {
        self.decoders.len()
    }
    pub fn is_empty(&self) -> bool {
        self.decoders.is_empty()
    }
}

/// The set of columns held by one chunk.
///
/// `BTreeMap`, not `HashMap`. Iteration order is part of the world hash, and a
/// `HashMap`'s order depends on a randomly-seeded hasher — which would make the
/// world hash differ between two runs of the *same binary*. That bug is trivial
/// to introduce and genuinely nasty to find; the fix is to make it impossible.
#[derive(Default)]
pub struct ColumnSet {
    columns: BTreeMap<ColumnId, Box<dyn Column>>,
}

impl Clone for ColumnSet {
    fn clone(&self) -> Self {
        ColumnSet {
            columns: self.columns.iter().map(|(k, v)| (*k, v.clone_box())).collect(),
        }
    }
}

impl ColumnSet {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn insert(&mut self, col: Box<dyn Column>) {
        self.columns.insert(col.id(), col);
    }

    pub fn get<T: Field + std::fmt::Debug>(&self, id: ColumnId) -> Option<&VecColumn<T>> {
        self.columns.get(&id)?.as_any().downcast_ref::<VecColumn<T>>()
    }

    pub fn get_mut<T: Field + std::fmt::Debug>(
        &mut self,
        id: ColumnId,
    ) -> Option<&mut VecColumn<T>> {
        self.columns.get_mut(&id)?.as_any_mut().downcast_mut::<VecColumn<T>>()
    }

    pub fn ids(&self) -> impl Iterator<Item = &ColumnId> {
        self.columns.keys()
    }

    pub fn len(&self) -> usize {
        self.columns.len()
    }
    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }

    pub fn encode(&self, w: &mut Writer) {
        w.u16(self.columns.len() as u16);
        for (id, col) in &self.columns {
            w.u16(id.0);
            col.encode(w);
        }
    }

    pub fn decode(r: &mut Reader, reg: &ColumnRegistry) -> Result<ColumnSet, CodecError> {
        let n = r.u16()?;
        let mut set = ColumnSet::new();
        for _ in 0..n {
            let id = ColumnId(r.u16()?);
            set.insert(reg.decode(id, r)?);
        }
        Ok(set)
    }
}

impl WorldHash for ColumnSet {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u16(self.columns.len() as u16);
        for col in self.columns.values() {
            col.hash_into(h);
        }
    }
}
