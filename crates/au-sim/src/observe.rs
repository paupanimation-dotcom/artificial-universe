//! Observation — a way to look at the world.
//!
//! # Why this exists, and why it is late
//!
//! Everything in this engine has been verified through tests. That has worked
//! well for *correctness* — the suite has caught a ledger netting opposing
//! flows, a vesicle dissolving itself, a reaction creating a carbon atom — but
//! it is close to useless for *discovery*, because a test can only check
//! something you already thought to ask.
//!
//! The artificial-life literature is unanimous on this point: finding
//! interesting behaviour requires playing with the thing, watching it, and
//! noticing what you did not predict. This project has had no way to do that,
//! and it has cost real time. Three separate debugging sessions ended with a
//! `eprintln!` in a hot loop and a `grep`, and each of them was solved the
//! moment actual numbers appeared:
//!
//!   * a world cycling every two ticks — visible instantly as a sawtooth in the
//!     protocell census;
//!   * bags whose contents never changed — a flat line;
//!   * a cell at 10¹⁸ K, then at absolute zero — a temperature trace leaving
//!     the axis;
//!   * compartments starved of substrate by their own medium — one series
//!     falling to zero while another stayed flat.
//!
//! Four time series would have shown all four in seconds. That is the whole
//! ambition here: not a renderer, not a game view, just **numbers over time,
//! plotted, in a file you can open**.
//!
//! # Read-only, and structurally so
//!
//! A recorder that perturbed the world would be worse than no recorder. Phase 1
//! shipped a bug where an event emitted during chunk eviction made *RAM
//! pressure* part of world identity — the world hashed differently depending on
//! what had been paged out. Observation must never do anything like that.
//!
//! So [`Recorder::observe`] takes `&World`, not `&mut World`. It cannot write to
//! a column, cannot emit an event, cannot touch the entity allocator, and cannot
//! change the hash, because the borrow checker will not let it. That is not a
//! convention to be careful about; it is a compile error.
//!
//! It also lives outside the scheduler. Nothing in the layer stack depends on
//! it, no system runs it, and a world observed a thousand times is byte-identical
//! to one never observed at all. There is a test for that, because "obviously
//! read-only" is exactly the kind of claim this project has learned to check.
//!
//! # What is recorded, and why these
//!
//! Not everything — a sample is cheap only if it stays cheap. The set below is
//! chosen because each entry corresponds to a mistake that has actually
//! happened:
//!
//!   * **bulk vs encapsulated populations, separately.** Since protocells,
//!     matter lives in two places, and every conservation confusion in Phase 5i
//!     came from conflating them. Seeing them apart is the point.
//!   * **temperature min/mean/max.** Two of the worst bugs announced themselves
//!     as absurd temperatures and nobody was looking.
//!   * **the protocell census and its lifetime counters.** Births, divisions and
//!     deaths distinguish a growing population from a churning one — which is
//!     precisely the distinction a division count alone hides.
//!   * **per-bag contents.** A bag that never changes is not alive, and no
//!     aggregate reveals that.

use crate::world::World;
use au_physics::{derive, Energy, MaterialId};

/// One bag at one instant.
#[derive(Clone, Debug, PartialEq)]
pub struct BagSample {
    pub id: u32,
    pub generation: u32,
    /// `EntityId::NONE.index` when this bag nucleated rather than being born.
    pub parent: u32,
    pub born_tick: u64,
    pub cell: u32,
    /// Total molecules held, of every kind.
    pub contents: i128,
}

/// The world at one instant, as far as an observer can see it.
#[derive(Clone, Debug, PartialEq)]
pub struct Sample {
    pub tick: u64,
    pub seconds: f64,

    /// Molecules of each species sitting in the medium.
    pub bulk: Vec<i128>,
    /// Molecules of each species held inside protocells.
    ///
    /// Kept apart from `bulk` deliberately. Their sum is what conservation
    /// checks; their *ratio* is what tells you whether compartments are doing
    /// anything at all.
    pub encapsulated: Vec<i128>,

    /// Kelvin, over cells that have mass.
    pub temp_min: f64,
    pub temp_mean: f64,
    pub temp_max: f64,

    pub protocells: usize,
    pub nucleated: u64,
    pub divided: u64,
    pub lysed: u64,

    /// Energy that has crossed the world boundary, in and out.
    pub energy_in: i128,
    pub energy_out: i128,
    /// Atoms across the boundary, per element (Phase 5h).
    pub atoms_in: Vec<i128>,
    pub atoms_out: Vec<i128>,

    /// A prefix of the living bags, capped so a large population cannot make
    /// observation the most expensive thing in the run.
    pub bags: Vec<BagSample>,
}

/// Samples a world at intervals, and never changes it.
#[derive(Clone, Debug)]
pub struct Recorder {
    every: u64,
    max_bags: usize,
    samples: Vec<Sample>,
}

impl Recorder {
    /// Sample every `every` ticks, keeping at most `max_bags` individual bags
    /// per sample.
    pub fn new(every: u64, max_bags: usize) -> Recorder {
        Recorder { every: every.max(1), max_bags, samples: Vec::new() }
    }

    pub fn samples(&self) -> &[Sample] {
        &self.samples
    }

    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Take a sample if this tick is due. Cheap to call every tick.
    ///
    /// `&World` — not `&mut`. See the module docs: this is the guarantee, and it
    /// is enforced by the compiler rather than by care.
    pub fn observe(&mut self, world: &World) {
        let tick = world.clock.tick().0;
        if tick % self.every != 0 && !self.samples.is_empty() {
            return;
        }
        if self.samples.last().map(|s| s.tick) == Some(tick) {
            return;
        }
        self.samples.push(sample(world, self.max_bags));
    }

    /// Sample unconditionally, whatever the tick.
    pub fn snapshot(&mut self, world: &World) {
        self.samples.push(sample(world, self.max_bags));
    }
}

/// Read the world once.
pub fn sample(world: &World, max_bags: usize) -> Sample {
    let n_species = world.chem_registry.len();
    let mut bulk = vec![0i128; n_species];
    let mut encapsulated = vec![0i128; n_species];

    let mut t_min = f64::INFINITY;
    let mut t_max = f64::NEG_INFINITY;
    let mut t_sum = 0.0f64;
    let mut t_count = 0usize;

    for coord in world.active.iter() {
        let Some(chunk) = world.chunks.get(*coord) else { continue };
        let cols = chunk.columns();

        for (s, slot) in bulk.iter_mut().enumerate() {
            if let Some(c) = cols.get::<i128>(crate::chem_columns::species_column(s)) {
                *slot += c.as_slice().iter().sum::<i128>();
            }
        }

        let (Some(mass), Some(energy), Some(matid)) = (
            cols.get::<i128>(crate::physics_columns::PHYS_MASS),
            cols.get::<i128>(crate::physics_columns::PHYS_ENERGY),
            cols.get::<u16>(crate::physics_columns::PHYS_MATERIAL),
        ) else {
            continue;
        };
        let (mass, energy, matid) = (mass.as_slice(), energy.as_slice(), matid.as_slice());
        for c in 0..mass.len().min(energy.len()).min(matid.len()) {
            let kg = mass[c] as f64 / au_physics::AG_PER_KG as f64;
            if kg <= 0.0 {
                continue;
            }
            let Some(m) = world.materials.get(MaterialId(matid[c])) else { continue };
            let st = derive(kg, Energy(energy[c]).as_joules(), m, 0.0, 1.0);
            let t = st.temperature;
            if t.is_finite() {
                t_min = t_min.min(t);
                t_max = t_max.max(t);
                t_sum += t;
                t_count += 1;
            }
        }
    }

    world.protocells.totals_into(&mut encapsulated);

    let bags: Vec<BagSample> = world
        .protocells
        .iter()
        .take(max_bags)
        .map(|p| BagSample {
            id: p.id.index,
            generation: p.id.generation,
            parent: p.parent.index,
            born_tick: p.born_tick,
            cell: p.cell,
            contents: p.total(),
        })
        .collect();

    Sample {
        tick: world.clock.tick().0,
        seconds: world.clock.now().as_secs_f64(),
        bulk,
        encapsulated,
        temp_min: if t_count > 0 { t_min } else { 0.0 },
        temp_mean: if t_count > 0 { t_sum / t_count as f64 } else { 0.0 },
        temp_max: if t_count > 0 { t_max } else { 0.0 },
        protocells: world.protocells.len(),
        nucleated: world.protocells.stat_nucleated,
        divided: world.protocells.stat_divided,
        lysed: world.protocells.stat_lysed,
        energy_in: world.ledger.energy_in.0,
        energy_out: world.ledger.energy_out.0,
        atoms_in: world.atoms.inflow.clone(),
        atoms_out: world.atoms.outflow.clone(),
        bags,
    }
}
