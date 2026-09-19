//! Where chemistry meets storage — the twin of `physics_columns`.
//!
//! `au-chem` knows nothing about chunks, and `au-data` knows nothing about
//! molecules. This file is the only place that knows both, and like its physics
//! counterpart it stays deliberately small: the moment it starts making chemical
//! decisions, the layering has sprung a leak.
//!
//! # The column-id range is permanent
//!
//! Physics owns `0x1xxx`; chemistry owns `0x2xxx`. A save file names its columns
//! by number, so renumbering one does not fail to load — it loads *wrong*,
//! silently reading one species' population as another's. These numbers are now
//! immortal.
//!
//! # What is stored per cell, and what is not
//!
//! Per cell, chemistry stores **one integer population per species** — a count of
//! molecules of that species in that cell, `i128`, exact, conserved. Nothing else.
//! In particular it does *not* store its own energy or temperature:
//!
//!   * **Energy is shared with physics.** Chemistry reads and writes the *same*
//!     `PHYS_ENERGY` column the thermodynamics solver conducts and radiates. That
//!     sharing is the whole point of Phase 3 — an exothermic reaction warms the
//!     cell that physics then carries heat away from. There is one energy account,
//!     not a chemical one beside a thermal one.
//!
//!   * **Temperature is derived, never stored.** Exactly as in physics: chemistry
//!     asks `au_physics::derive` what temperature the shared energy implies, uses
//!     it to set Arrhenius rates, and writes back only energy. Temperature has no
//!     column here for the same reason it has none there — storing it would let it
//!     drift out of agreement with the energy it is supposed to be a function of.
//!
//! # Why the species count is fixed at boot (and why that is not cheating)
//!
//! The columns below are numbered `CHEM_SPECIES_BASE + i`, one per species, and
//! the *number* of species a world has is fixed when the world boots — declared in
//! config, the same way the periodic table and the material table are declared.
//!
//! This is the honest scope of Phase 3b. Open-ended molecular discovery — where a
//! reaction invents a graph no config ever named, and the cell grows a new
//! population slot for it on the fly — is real and is coming, but it is a Phase 5+
//! problem (it needs the species registry to become serialized state and the
//! per-cell storage to become sparse). Fixing the *vocabulary* of molecules at
//! boot does not fix the *outcome*: which reactions fire, what equilibria they
//! reach, how much heat they release, whether a cell ignites — all of that is
//! emergent from rules interacting. Only the list of molecules that *can* exist is
//! declared, precisely as "which elements exist" is declared. Nothing about a wolf
//! is written down; only that carbon is available.

use au_data::column::{ColumnId, ColumnRegistry, ColumnSet, VecColumn};

/// The first chemistry column id. Species `i` lives at `CHEM_SPECIES_BASE + i`.
///
/// `0x2000` begins the chemistry range. With this base, a world may declare up to
/// `0x1000` = 4096 species before colliding with anything — far more than Phase 3b
/// needs, and the ceiling moves when genetics (`0x3xxx`) is still distant.
pub const CHEM_SPECIES_BASE: u16 = 0x2000;

/// The maximum number of species a Phase 3b world may declare. A guard, not a
/// prophecy: it keeps the chemistry range from ever reaching into a future layer's.
pub const CHEM_MAX_SPECIES: usize = 0x0F00; // 3840

/// The column id for species index `i`.
#[inline]
pub fn species_column(i: usize) -> ColumnId {
    debug_assert!(i < CHEM_MAX_SPECIES, "species index {} exceeds the chemistry range", i);
    ColumnId(CHEM_SPECIES_BASE + i as u16)
}

/// Register `n_species` population columns. Called once at world construction,
/// after the species count is known from config, so the snapshot codec can
/// round-trip them through the same generic path physics columns use.
pub fn register(reg: &mut ColumnRegistry, n_species: usize) {
    assert!(
        n_species <= CHEM_MAX_SPECIES,
        "a world declared {} species; the chemistry column range holds {}",
        n_species,
        CHEM_MAX_SPECIES
    );
    for i in 0..n_species {
        reg.register::<i128>(species_column(i));
    }
}

/// Does this chunk carry chemistry — at least the first species column?
///
/// Most chunks will not. Chemistry, like physics, runs on a thin, deliberately
/// small active set, not on the vacuum that is most of a world.
pub fn has_chemistry(cols: &ColumnSet, n_species: usize) -> bool {
    n_species > 0 && cols.get::<i128>(species_column(0)).is_some()
}

/// Install empty (all-zero) population columns for `n_species` species.
///
/// Called by generators, never by systems — a system that could *create* a
/// population column where there was none is a system that can conjure matter, and
/// conjuring is the thing the Golden Rule exists to forbid. A freshly installed
/// chemistry chunk is empty: no molecules, which is the correct starting point for
/// a cell whose contents must arrive by cause, not by fiat.
pub fn install(cols: &mut ColumnSet, cells: usize, n_species: usize) {
    for i in 0..n_species {
        cols.insert(Box::new(VecColumn::<i128>::with_data(species_column(i), vec![0; cells])));
    }
}
