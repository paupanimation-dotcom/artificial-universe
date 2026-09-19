//! Stable identity.
//!
//! Organisms die. Their slots get reused. If identity were just an index, a
//! family tree would eventually point at a stranger — and family trees are a
//! feature we promised the player (VISION, Player Experience). So identity is
//! index + generation: reusing a slot bumps the generation, and a stale handle
//! is detectably stale rather than silently wrong.

use crate::hash::{Hasher, WorldHash};

/// A handle to a simulation entity. 8 bytes, `Copy`, dense.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub struct EntityId {
    pub index: u32,
    pub generation: u32,
}

impl EntityId {
    pub const NONE: EntityId = EntityId { index: u32::MAX, generation: u32::MAX };

    pub fn is_none(self) -> bool {
        self == EntityId::NONE
    }
}

impl WorldHash for EntityId {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u32(self.index);
        h.write_u32(self.generation);
    }
}

/// Allocates entity slots and recycles dead ones.
///
/// Deterministic by construction: the free list is LIFO and is part of the
/// snapshot, so a reloaded world hands out the same ids as a continuous one.
#[derive(Clone, Debug, Default)]
pub struct EntityAllocator {
    generations: Vec<u32>,
    free: Vec<u32>,
}

impl EntityAllocator {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn alloc(&mut self) -> EntityId {
        if let Some(index) = self.free.pop() {
            EntityId { index, generation: self.generations[index as usize] }
        } else {
            let index = self.generations.len() as u32;
            self.generations.push(0);
            EntityId { index, generation: 0 }
        }
    }

    pub fn free(&mut self, id: EntityId) -> bool {
        if !self.is_live(id) {
            return false;
        }
        self.generations[id.index as usize] = self.generations[id.index as usize].wrapping_add(1);
        self.free.push(id.index);
        true
    }

    pub fn is_live(&self, id: EntityId) -> bool {
        (id.index as usize) < self.generations.len()
            && self.generations[id.index as usize] == id.generation
            && !self.free.contains(&id.index)
    }

    pub fn capacity(&self) -> usize {
        self.generations.len()
    }

    pub fn live_count(&self) -> usize {
        self.generations.len() - self.free.len()
    }

    pub fn raw(&self) -> (&[u32], &[u32]) {
        (&self.generations, &self.free)
    }

    pub fn restore(generations: Vec<u32>, free: Vec<u32>) -> Self {
        EntityAllocator { generations, free }
    }
}

impl WorldHash for EntityAllocator {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u64(self.generations.len() as u64);
        for &g in &self.generations {
            h.write_u32(g);
        }
        h.write_u64(self.free.len() as u64);
        for &f in &self.free {
            h.write_u32(f);
        }
    }
}
