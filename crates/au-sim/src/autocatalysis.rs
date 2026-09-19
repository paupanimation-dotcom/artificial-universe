//! Looking at a world's chemistry and asking whether it has closed on itself.
//!
//! # Read-only, on purpose
//!
//! Everything in this module is a *query*. It reads the registry and the
//! populations, works on a **clone** of the registry so that even the act of
//! enumerating cannot intern a species into the real world, and returns a
//! report. It emits no events, writes no columns, and touches no clock. Delete
//! the module and every world runs bit-identically.
//!
//! That restraint is deliberate and it is not merely tidiness. An event in the
//! log is world state — it is hashed, snapshotted, and resumed — so making "a
//! RAF appeared" an event would put an *observation* inside the thing observed.
//! Worse, it would not survive resume: the emit-once rule depends on knowing
//! that the previous tick had no RAF, and that fact is not in any snapshot, so a
//! resumed world would re-announce a discovery it had already made and diverge
//! from one that never stopped. Phase 5a's discovery events are safe precisely
//! because they are tied to registry growth, which *is* serialized.
//!
//! When a timeline system arrives, "the chemistry first closed on itself" is an
//! obvious thing for it to record — but it will have to be recorded from state
//! that a snapshot carries. Until then this stays a lens.
//!
//! # The food set, and the modelling choice in it
//!
//! A RAF is defined relative to a **food set**: the species the environment
//! supplies for free. This module takes the food to be the **monatomic
//! species** — the bare elements. That is the geochemically meaningful reading:
//! volcanism, weathering and infall supply atoms, and the question worth asking
//! is whether a network can sustain itself *given elemental feedstock*.
//!
//! The limitation is worth stating plainly: a sealed demo cell does not
//! actually replenish its atoms, so a RAF found there answers a slightly
//! idealised question — "could this network sustain itself if fed?" — rather
//! than "is it sustaining itself right now". Sustained-in-fact requires flux
//! boundary conditions the chemistry does not yet have.

use au_chem::{
    catalysed_variants, enumerate_reactions, maximal_raf, self_producing, RafReport, Reaction,
    SpeciesId,
};

use crate::chem_columns::{has_chemistry, species_column};
use crate::chemistry::Chemistry;
use crate::world::World;

/// What the lens saw, with enough context to print it.
pub struct Survey {
    /// The bare (uncatalysed) network the survey was run on.
    pub reactions: Vec<Reaction>,
    /// The catalysis relation: (catalyst, index into `reactions`).
    pub relation: Vec<(SpeciesId, usize)>,
    /// The food set used — the monatomic species.
    pub food: Vec<SpeciesId>,
    /// Every species with a nonzero population somewhere, plus the seeds.
    pub present: Vec<SpeciesId>,
    /// The maximal RAF, if the network has one.
    pub raf: Option<RafReport>,
    /// Species inside the RAF that catalyse a reaction producing them — the
    /// ones that are, in the plainest sense, making more of themselves.
    pub self_producing: Vec<SpeciesId>,
}

impl Survey {
    /// How many distinct species act as catalysts anywhere in the network.
    pub fn catalyst_count(&self) -> usize {
        let mut v: Vec<u32> = self.relation.iter().map(|(c, _)| c.0).collect();
        v.sort_unstable();
        v.dedup();
        v.len()
    }
}

/// Survey a world's chemistry for autocatalytic sets.
///
/// Returns `None` when the world has no open chemistry to survey — a declared
/// Phase 3 reaction list has no structure to derive catalysis from, so the
/// question is not applicable rather than answered in the negative.
pub fn survey(world: &World) -> Option<Survey> {
    survey_with(world, None)
}

/// Survey a world while *overriding* its catalysis rules.
///
/// The world is not modified and does not run — this asks a counterfactual
/// question about the network it has already discovered: what would close on
/// itself if a grip had to be this polar to count? Sweeping the override is how
/// the density of catalysis is walked down over a single fixed network, which is
/// the only way to see the collapse cleanly. Comparing separately-run worlds
/// would confound the change in catalysis with the change in what each world
/// happened to discover.
pub fn survey_with(world: &World, over: Option<au_chem::CatalysisRules>) -> Option<Survey> {
    let chem = match Chemistry::from_kv(world.config.as_map()) {
        Ok(Some(c)) => c,
        _ => return None,
    };
    let rules = chem.open?;
    let pool = chem.column_count();

    // Present = seeds ∪ everything with a population, exactly as the system
    // computes it, so the survey describes the network the world is actually
    // running rather than a hypothetical one.
    let mut present: std::collections::BTreeSet<u32> =
        (0..chem.species_count() as u32).collect();
    for coord in world.active.iter() {
        let Some(chunk) = world.chunks.get(*coord) else { continue };
        if !has_chemistry(chunk.columns(), pool) {
            continue;
        }
        for sp in 0..pool {
            let col = chunk.columns().get::<i128>(species_column(sp)).unwrap();
            if col.as_slice().iter().any(|&v| v != 0) {
                present.insert(sp as u32);
            }
        }
    }
    let ids: Vec<SpeciesId> = present.iter().map(|&i| SpeciesId(i)).collect();

    // The clone is the whole point: enumeration interns, and the instrument may
    // not change what it measures.
    let mut reg = world.chem_registry.clone();
    let reactions: Vec<Reaction> =
        enumerate_reactions(&ids, &mut reg, &chem.table, &chem.bond_model, &rules)
            .into_iter()
            .filter(|r| {
                r.products
                    .iter()
                    .chain(r.reactants.iter())
                    .all(|t| (t.species.0 as usize) < pool)
            })
            .collect();

    let relation: Vec<(SpeciesId, usize)> = match over.or(chem.catalysis) {
        Some(crules) => catalysed_variants(
            &reactions,
            &ids,
            &reg,
            &chem.table,
            &chem.bond_model,
            &crules,
            &rules,
        )
        .into_iter()
        .map(|v| (v.catalyst, v.parent))
        .collect(),
        None => Vec::new(),
    };

    // Food: the bare elements the environment is taken to supply.
    let food: Vec<SpeciesId> = (0..pool.min(reg.len()))
        .map(|i| SpeciesId(i as u32))
        .filter(|&s| reg.get(s).map(|m| m.atom_count() == 1).unwrap_or(false))
        .collect();

    let raf = maximal_raf(&reactions, &relation, &food);
    let sp = match &raf {
        Some(r) => self_producing(&reactions, &relation, r),
        None => Vec::new(),
    };

    Some(Survey { reactions, relation, food, present: ids, raf, self_producing: sp })
}
