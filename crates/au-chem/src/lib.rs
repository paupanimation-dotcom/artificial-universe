//! # au-chem — chemistry as a computable abstraction
//!
//! ## The honest framing
//!
//! Real chemistry is not computable at the scale this project runs at. Solving the
//! Schrödinger equation for molecules costs roughly O(N⁷) and runs on
//! supercomputers for *single* reactions; this simulation will have millions of
//! cells. There is no version of this crate that simulates literal electrons and
//! reaches biology. Every serious artificial-life project that promised "emergent
//! real chemistry" either never shipped it or quietly faked it.
//!
//! So this crate does not simulate chemistry. It provides an **abstraction of
//! chemistry chosen to preserve the causal structure that lets chemistry produce
//! life**, while staying computable. That is the actual frontier of Phase 3, and
//! the whole tower above it rests on the abstraction being right: detailed enough
//! that novelty and self-organisation can emerge, cheap enough to run a planet.
//!
//! ## The three commitments
//!
//! 1. **Elements are exactly conserved integer counts** ([`element`]). Atoms are
//!    rearranged, never created — because a chemistry that can leak a carbon atom
//!    is a chemistry something will evolve to farm.
//!
//! 2. **Molecules are graphs, discovered not listed** ([`molecule`]). Water is
//!    "two H bonded to one O", a graph the engine builds and fingerprints, not a
//!    row in a table. Genuine novelty is possible because a new molecule is just a
//!    new graph.
//!
//! 3. **Reactions are physics-gated rules** ([`reaction`], arriving next). A
//!    reaction fires by mass-action kinetics at an Arrhenius rate, and its energy
//!    of reaction flows into the Phase 2 energy field — so an exothermic reaction
//!    heats the cell it happens in, and chemistry and thermodynamics become one
//!    conserved loop rather than two systems bolted together.
//!
//! This crate depends on `au-physics` for exactly that reason: chemistry's energy
//! *is* physics' energy, drawn from the same conserved books.

pub mod amphiphile;
pub mod assembly;
pub mod catalysis;
pub mod element;
pub mod kinetics;
pub mod molecule;
pub mod network;
pub mod raf;
pub mod reaction;
pub mod reservoir;
pub mod vesicle;
pub mod reactor;

pub use amphiphile::{amphiphilicity, atom_polarity, best_cut, Cut};
pub use assembly::{
    assembled_count, coverage, critical_concentration, enclosing_area, membrane_molecules,
    permeability, score_table, MembraneRules,
};
pub use catalysis::{catalysed_variants, catalysts_for, catalytic_reduction, CatalysisRules, Catalysed};
pub use element::{AtomCount, Element, PeriodicTable, Z};
pub use kinetics::{
    bimolecular_prefactor, prefactor_at, termolecular_prefactor, unimolecular_prefactor,
    KineticRules,
};
pub use molecule::{Atom, Bond, BondOrder, CanonicalForm, Molecule};
pub use raf::{closure, maximal_raf, self_producing, RafReport};
pub use reservoir::{formula_table, port_delta, AtomLedger, Port};
pub use vesicle::{
    composition_distance, fate, partition, reduced_volume, shape, Fate, Shape,
    FISSION_REDUCED_VOLUME,
};
pub use reaction::{BondEnergyModel, Reaction, SpeciesId, SpeciesRegistry, Term, R_GAS};
pub use network::{
    structural_prefactor,
    enumerate_reactions, formula_string, join, split, with_order, NetworkRules, AVOGADRO,
};
pub use reactor::{react, with_derived_enthalpy, CellChemistry, ReactingCell, ReactorReport};
