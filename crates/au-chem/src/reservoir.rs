//! Open boundaries for matter — the ports through which atoms enter and leave.
//!
//! # Why a sealed world cannot hold life
//!
//! Energy has been open since Phase 2: heat enters at the floor, radiates out at
//! the top, and everything interesting lives in the gap. Matter has been sealed
//! the entire time. Every demo this engine has ever run — the reactor, the vent,
//! the RAF survey — is a box with a fixed number of atoms in it.
//!
//! A box has one destiny. It runs down to equilibrium and stays there. The
//! autocatalytic set Phase 5c found is real, and in a box it is also doomed: it
//! eats its food, makes its products, and stops. That is why the Phase 5c
//! changelog was careful to say the survey answers *"could this network sustain
//! itself if it were fed?"* rather than *"is it sustaining itself?"*. This module
//! is what turns the first question into the second.
//!
//! Life is a dissipative structure. Dissipative structures require throughput —
//! not merely a gradient in energy, but a *flow of matter* through the system.
//! You cannot have selection in a box that is dying, because every lineage in it
//! has the same fate.
//!
//! # Dirichlet, not injection
//!
//! The obvious implementation is a source term: add N atoms per second. It is
//! also the wrong one, twice over. It needs a fractional accumulator to avoid a
//! deadband at low rates, and — worse — it *declares the flux*, when flux is
//! exactly the thing that ought to be derived.
//!
//! So a port is a **Dirichlet condition on population**, the direct twin of the
//! `bottom_temp` / `top_temp` conditions the fluid solver has used since Phase 2b:
//!
//!   * a [`Port::Source`] is a window onto a reservoir held at a fixed
//!     composition — a vent in contact with something enormous;
//!   * a [`Port::Sink`] is a window onto an infinitely dilute exterior — an open
//!     port to the ocean.
//!
//! Nothing declares how fast matter moves. The port declares a *concentration*,
//! the gradient does the rest, and the flux that results is an outcome of the
//! geometry and the diffusivity. Every operation is an exact integer clamp, so
//! there is no rounding, no residue to carry, and no new state to snapshot.
//!
//! # A sink may not have taste
//!
//! [`Port::Sink`] draws **every** species to zero, and this is the single most
//! important constraint in the module. A sink that removed some molecules and
//! spared others would be a fitness function written into the boundary
//! condition — the exact thing this project exists not to do. The outflow is
//! blind. Whether anything survives it must be decided by what the chemistry
//! does, not by what the drain prefers.
//!
//! What makes that interesting is that the engine already has a way for a cell
//! to resist a blind drain: a membrane. Since Phase 5g, assembled surfactant
//! lowers the permeability of a cell's faces, and permeability multiplies the
//! diffusive exchange. So a cell whose own chemistry builds surfactant loses its
//! contents to the port more slowly than one that does not — and it does so
//! without anything in the code granting membranes an advantage. The advantage
//! is what a barrier *is*, in the presence of a flow. That is selection arriving
//! as a consequence of physics rather than as a parameter.
//!
//! # What the books must say
//!
//! `au-physics` already books every joule crossing a boundary, because a world
//! that can create energy will eventually evolve something that lives on the
//! leak. Matter now needs the same treatment, and per element rather than in
//! bulk: a bug that turned carbon into oxygen would conserve mass to the
//! attogram and still be a catastrophe. Hence [`AtomLedger`], and the identity
//! it exists to make checkable:
//!
//! ```text
//!     atoms_of_element_e(now) - atoms_of_element_e(start) == in[e] - out[e]
//! ```
//!
//! exactly, for every element, on every tick. Conservation does not weaken when
//! the world opens; it becomes a statement about bookkeeping instead of a
//! statement about a constant.

use crate::element::PeriodicTable;
use crate::molecule::Molecule;
use crate::reaction::{SpeciesId, SpeciesRegistry};
use au_core::hash::{Hasher, WorldHash};

// ═══════════════════════════════════════════════════════════════════════════
//  The books
// ═══════════════════════════════════════════════════════════════════════════

/// Atoms that have crossed the boundary of the world, per element.
///
/// Indexed by *element index* in the [`PeriodicTable`] — the same index
/// [`Molecule::formula`] returns — not by atomic number, so the vectors stay
/// dense however sparse the chosen elements are.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct AtomLedger {
    pub inflow: Vec<i128>,
    pub outflow: Vec<i128>,
}

impl AtomLedger {
    pub fn new(n_elements: usize) -> AtomLedger {
        AtomLedger { inflow: vec![0; n_elements], outflow: vec![0; n_elements] }
    }

    pub fn len(&self) -> usize {
        self.inflow.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inflow.is_empty()
    }

    /// Net atoms of element `e` the world should have gained since the start.
    pub fn net(&self, e: usize) -> i128 {
        self.inflow.get(e).copied().unwrap_or(0) - self.outflow.get(e).copied().unwrap_or(0)
    }

    /// Has anything crossed at all? A closed world answers `false` forever, and
    /// that is what lets the pre-5h invariant survive verbatim: where nothing
    /// crossed, `in - out == 0`, and "totals are constant" is the special case
    /// of the general identity rather than a rule that got replaced.
    pub fn any_flow(&self) -> bool {
        self.inflow.iter().chain(self.outflow.iter()).any(|&v| v != 0)
    }

    /// Book gross flows: molecules gained and molecules lost, separately.
    ///
    /// **Separately is the whole point, and it took a bug to make that obvious.**
    /// The first version of this accumulated one signed delta per species and
    /// booked that. It looks equivalent and is not: at steady state a source
    /// feeding 2,490 molecules and a sink draining 2,490 sum to zero, so the
    /// books recorded a sealed world while matter poured through it at full rate.
    /// The net identity still held — `here(now) - here(start) == in - out` was
    /// true, because both sides were zero — which is exactly what made it
    /// dangerous: the invariant everything else is checked against cannot detect
    /// this failure, and the only visible symptom was a throughput of nothing in
    /// a world that was obviously flowing.
    ///
    /// A ledger that nets opposing flows cannot answer the question ledgers exist
    /// for. Gross in, gross out, always.
    pub fn book_gross(&mut self, gained: &[i128], lost: &[i128], formula: &[Vec<i128>]) {
        for (s, f) in formula.iter().enumerate() {
            let g = gained.get(s).copied().unwrap_or(0);
            let l = lost.get(s).copied().unwrap_or(0);
            if g == 0 && l == 0 {
                continue;
            }
            for (e, &count) in f.iter().enumerate() {
                if count == 0 {
                    continue;
                }
                self.inflow[e] += count * g;
                self.outflow[e] += count * l;
            }
        }
    }

    /// Book a change in population as atoms crossing the boundary.
    ///
    /// `delta[s]` is the change a port made to species `s`: positive means the
    /// world gained molecules (inflow), negative means it lost them (outflow).
    /// `formula[s][e]` is the atom count of element `e` in species `s`,
    /// precomputed by [`formula_table`] because doing it per tick would mean
    /// walking every molecular graph in the registry every tick.
    pub fn book(&mut self, delta: &[i128], formula: &[Vec<i128>]) {
        for (s, &d) in delta.iter().enumerate() {
            if d == 0 {
                continue;
            }
            let Some(f) = formula.get(s) else { continue };
            for (e, &count) in f.iter().enumerate() {
                if count == 0 {
                    continue;
                }
                let atoms = count * d.abs();
                if d > 0 {
                    self.inflow[e] += atoms;
                } else {
                    self.outflow[e] += atoms;
                }
            }
        }
    }
}

impl WorldHash for AtomLedger {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_u32(self.inflow.len() as u32);
        for v in self.inflow.iter().chain(self.outflow.iter()) {
            h.write_bytes(&v.to_le_bytes());
        }
    }
}

/// Atom counts per element for every species in the registry, in species order.
///
/// Species discovered later append to the registry, so this table is rebuilt
/// whenever the registry grows — the same cadence as the reaction cache.
pub fn formula_table(reg: &SpeciesRegistry, table: &PeriodicTable) -> Vec<Vec<i128>> {
    (0..reg.len())
        .map(|i| {
            reg.get(SpeciesId(i as u32))
                .map(|m: &Molecule| m.formula(table))
                .unwrap_or_else(|| vec![0; table.len()])
        })
        .collect()
}

// ═══════════════════════════════════════════════════════════════════════════
//  Ports
// ═══════════════════════════════════════════════════════════════════════════

/// A boundary condition on the species populations of one cell.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Port {
    /// A window onto a reservoir of fixed composition.
    ///
    /// Each listed species is topped **up** to its stated population and never
    /// down. A source that also removed would be a sink wearing a source's
    /// clothes, and the two would become impossible to tell apart in the books.
    /// Species not listed are untouched — a vent supplies what it supplies and
    /// is indifferent to everything else in the water.
    Source { hold: Vec<(SpeciesId, i128)> },

    /// A window onto an infinitely dilute exterior: every species drawn to zero.
    ///
    /// Blind by construction. See the module docs for why that matters more than
    /// anything else here.
    Sink,
}

impl Port {
    /// Does this source supply nothing but bare elements?
    ///
    /// Not enforced by the type, because a real hydrothermal vent emits H₂, CH₄
    /// and H₂S rather than a mist of loose atoms, and forbidding that would be
    /// dogma rather than physics. But it is worth being able to *ask*, because
    /// there is one genuine trap here: **a boundary that supplies a member of
    /// the set under study is not an experiment, it is an assumption.** Feed a
    /// world the very molecule whose emergence you are trying to observe and you
    /// have proven nothing. The chemostat demo asserts this of its own feed.
    pub fn is_elemental(&self, reg: &SpeciesRegistry) -> bool {
        match self {
            Port::Sink => true,
            Port::Source { hold } => hold
                .iter()
                .all(|(s, _)| reg.get(*s).map(|m| m.atom_count() == 1).unwrap_or(false)),
        }
    }
}

/// What a port would do to a cell, as a change per species.
///
/// Nothing is written here: the caller owns the storage and the layout (the
/// simulation holds populations species-major, which is the wrong shape to hand
/// over as a slice). `pop_at(s)` reads the cell's current population of species
/// `s`; `delta` is resized to `n_species` and filled with the change required.
///
/// Pure, total, and exact — every value is an integer clamp, so applying the
/// same port twice in a row is a no-op the second time.
pub fn port_delta(
    port: &Port,
    n_species: usize,
    pop_at: impl Fn(usize) -> i128,
    delta: &mut Vec<i128>,
) {
    delta.clear();
    delta.resize(n_species, 0);
    match port {
        Port::Source { hold } => {
            for &(s, target) in hold {
                let i = s.0 as usize;
                if i >= n_species {
                    continue;
                }
                let have = pop_at(i);
                if have < target {
                    delta[i] = target - have;
                }
            }
        }
        Port::Sink => {
            for (i, d) in delta.iter_mut().enumerate().take(n_species) {
                let have = pop_at(i);
                if have != 0 {
                    *d = -have;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(v: &[i128]) -> impl Fn(usize) -> i128 + '_ {
        move |i| v.get(i).copied().unwrap_or(0)
    }

    #[test]
    fn a_source_tops_up_and_never_removes() {
        let pops = [500i128, 9_000, 0];
        let port = Port::Source {
            hold: vec![(SpeciesId(0), 1_000), (SpeciesId(1), 1_000), (SpeciesId(2), 1_000)],
        };
        let mut d = Vec::new();
        port_delta(&port, 3, cell(&pops), &mut d);
        assert_eq!(d, vec![500, 0, 1_000], "an over-full cell must not be drained");
    }

    #[test]
    fn a_sink_takes_everything_without_preference() {
        let pops = [7i128, 0, 123_456];
        let mut d = Vec::new();
        port_delta(&Port::Sink, 3, cell(&pops), &mut d);
        assert_eq!(d, vec![-7, 0, -123_456]);
    }

    #[test]
    fn applying_a_port_twice_changes_nothing_the_second_time() {
        let mut pops = vec![0i128, 40];
        let port = Port::Source { hold: vec![(SpeciesId(0), 100), (SpeciesId(1), 100)] };
        let mut d = Vec::new();
        port_delta(&port, 2, cell(&pops), &mut d);
        for (p, dd) in pops.iter_mut().zip(&d) {
            *p += dd;
        }
        port_delta(&port, 2, cell(&pops), &mut d);
        assert_eq!(d, vec![0, 0]);
    }

    #[test]
    fn the_ledger_books_outflow_by_element() {
        // Species 0 is a bare atom of element 0; species 1 is E0-E1-E1.
        let formula = vec![vec![1, 0], vec![1, 2]];
        let mut l = AtomLedger::new(2);
        l.book(&[-3, -5], &formula);
        assert_eq!(l.outflow, vec![8, 10]);
        assert_eq!(l.inflow, vec![0, 0]);
        assert_eq!(l.net(0), -8);
    }

    /// The bug that made `book_gross` exist. Two ports, opposite directions,
    /// equal size — the state of a world at steady state — must be recorded as
    /// throughput, not as nothing happening.
    #[test]
    fn opposing_ports_are_booked_gross_not_net() {
        let formula = vec![vec![1, 0]];
        let mut l = AtomLedger::new(2);
        l.book_gross(&[2_490], &[2_490], &formula);
        assert_eq!(l.inflow[0], 2_490);
        assert_eq!(l.outflow[0], 2_490);
        assert_eq!(l.net(0), 0, "the net is still zero — that part was never wrong");
        assert!(l.any_flow(), "a world at steady state is not a sealed one");
    }

    #[test]
    fn a_closed_world_reports_no_flow() {
        let l = AtomLedger::new(3);
        assert!(!l.any_flow());
    }
}
