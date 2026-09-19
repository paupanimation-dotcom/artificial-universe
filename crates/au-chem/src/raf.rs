//! Finding the sets of reactions that collectively make themselves.
//!
//! # An instrument, not a mechanism
//!
//! This module changes nothing. It runs no chemistry, moves no atoms, alters no
//! rate. It is a *lens* pointed at the reaction network a world has already
//! discovered, and it answers one question: is there a set of reactions that
//! collectively catalyses itself, using only what the environment freely
//! supplies?
//!
//! That distinction is the whole reason this is allowed to exist. The project's
//! rule is that nothing important may be hand-placed — so the engine may not
//! *create* a replicator. It may, however, **notice** one, exactly as the
//! timeline system notices an extinction without causing it. Detection is
//! observation. If this module is deleted, every world runs identically.
//!
//! # RAF: what "collectively makes itself" means precisely
//!
//! The definition is Hordijk & Steel's (2004), building on Kauffman's
//! autocatalytic sets (1971, 1986). Given a set of reactions, a catalysis
//! relation, and a **food set** F of species the environment supplies for free,
//! a subset R′ of the reactions is a **RAF** — Reflexively Autocatalytic and
//! F-generated — when both of these hold:
//!
//!   * **RA (reflexively autocatalytic).** Every reaction in R′ is catalysed by
//!     at least one species that R′ itself can produce from F. No outside help:
//!     the set supplies its own catalysts.
//!   * **F-generated.** Every reactant of every reaction in R′ is either in F or
//!     is built by other reactions in R′. No outside help: the set supplies its
//!     own substrates, given food.
//!
//! Together: a closed loop of production in which everything needed to run the
//! loop is made by the loop, out of what the world hands over freely. That is
//! not life, and this module never says it is — a RAF has no membrane, no
//! heredity, no individuality, and cannot evolve. It is the *network-level
//! precondition* for those things: metabolism before there is anything for the
//! metabolism to belong to.
//!
//! # The algorithm, and why it terminates
//!
//! Finding the maximal RAF looks like it should be exponential — there are 2^|R|
//! subsets — and it is not, which is the elegant part. Start with every reaction
//! and prune:
//!
//! ```text
//!   R′ := all reactions
//!   loop:
//!       W  := closure of F under R′       (everything R′ can build from food)
//!       R″ := { r ∈ R′ : r's reactants ⊆ W, and some catalyst of r is in W }
//!       if R″ = R′ : stop
//!       R′ := R″
//! ```
//!
//! Each round removes at least one reaction or halts, so there are at most |R|
//! rounds, each of polynomial cost. The result is *the* maximal RAF: pruning can
//! only ever remove a reaction that no RAF could contain, so nothing that
//! belongs is lost, and the answer is unique — the union of all RAFs is itself a
//! RAF. An empty result is a real answer: this network has no self-sustaining
//! core.
//!
//! Everything here is sets and integers in sorted order. No floats, no RNG, no
//! iteration over a hash map. The same network gives the same RAF on every
//! machine, forever, which is the only way an observation can be part of a
//! world's history.

use crate::reaction::{Reaction, SpeciesId};
use std::collections::BTreeSet;

/// What the lens saw.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RafReport {
    /// Indices into the reaction slice that survived pruning, ascending.
    pub reactions: Vec<usize>,
    /// The closure: every species the surviving set can build from the food,
    /// ascending. Includes the food itself.
    pub closure: Vec<SpeciesId>,
    /// How many pruning rounds it took. Diagnostic only.
    pub rounds: u32,
}

/// Everything reachable from `food` using only the reactions marked active:
/// the smallest set containing the food and closed under the active reactions.
pub fn closure(reactions: &[Reaction], active: &[bool], food: &BTreeSet<SpeciesId>) -> BTreeSet<SpeciesId> {
    let mut w = food.clone();
    loop {
        let mut grew = false;
        for (i, r) in reactions.iter().enumerate() {
            if !active[i] {
                continue;
            }
            if r.reactants.iter().all(|t| w.contains(&t.species)) {
                for t in &r.products {
                    if w.insert(t.species) {
                        grew = true;
                    }
                }
            }
        }
        if !grew {
            return w;
        }
    }
}

/// The maximal RAF of a catalytic reaction system, or `None` if there is none.
///
/// `catalysis` is a list of (catalyst species, reaction index) pairs — the
/// relation, kept separate from the reactions themselves. This separation is
/// deliberate and it matters: the *kinetics* run catalysed reactions with the
/// catalyst written on both sides (see [`crate::catalysis`]), but a RAF must be
/// computed on the bare network plus the relation. Fed the catalysed variants
/// instead, the "is it catalysed by something in the closure" test would be
/// trivially satisfied by the catalyst sitting among the reactants, and every
/// answer would be yes.
pub fn maximal_raf(
    reactions: &[Reaction],
    catalysis: &[(SpeciesId, usize)],
    food: &[SpeciesId],
) -> Option<RafReport> {
    let n = reactions.len();
    if n == 0 {
        return None;
    }
    let food: BTreeSet<SpeciesId> = food.iter().copied().collect();

    // Catalysts of each reaction, gathered once.
    let mut cats: Vec<BTreeSet<SpeciesId>> = vec![BTreeSet::new(); n];
    for &(c, i) in catalysis {
        if i < n {
            cats[i].insert(c);
        }
    }

    let mut active = vec![true; n];
    let mut rounds = 0u32;
    loop {
        rounds += 1;
        let w = closure(reactions, &active, &food);
        let mut changed = false;
        for i in 0..n {
            if !active[i] {
                continue;
            }
            let substrates_ok = reactions[i].reactants.iter().all(|t| w.contains(&t.species));
            let catalysed = cats[i].iter().any(|c| w.contains(c));
            if !(substrates_ok && catalysed) {
                active[i] = false;
                changed = true;
            }
        }
        if !changed {
            let surviving: Vec<usize> = (0..n).filter(|&i| active[i]).collect();
            if surviving.is_empty() {
                return None;
            }
            // Report the closure of the surviving set, which may be smaller than
            // the last one computed if the final round pruned nothing but an
            // earlier one did.
            let w = closure(reactions, &active, &food);
            return Some(RafReport {
                reactions: surviving,
                closure: w.into_iter().collect(),
                rounds,
            });
        }
    }
}

/// The species in a RAF that are *stoichiometrically* autocatalytic under
/// catalysis: species which catalyse a reaction that produces them.
///
/// This is the sharpest thing the lens can report. Such a species is, in the
/// plainest possible sense, making more of itself: its presence speeds the very
/// reaction whose output it is. Phase 5a's network could not contain one — the
/// derived moves have no room for a species on both sides — so finding one is
/// always a fact about catalysis having opened the door.
pub fn self_producing(
    reactions: &[Reaction],
    catalysis: &[(SpeciesId, usize)],
    within: &RafReport,
) -> Vec<SpeciesId> {
    let inside: BTreeSet<usize> = within.reactions.iter().copied().collect();
    let mut out = BTreeSet::new();
    for &(c, i) in catalysis {
        if !inside.contains(&i) || i >= reactions.len() {
            continue;
        }
        if reactions[i].products.iter().any(|t| t.species == c) {
            out.insert(c);
        }
    }
    out.into_iter().collect()
}
