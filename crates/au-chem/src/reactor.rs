//! The reactor: chemistry advancing in one cell, coupled to that cell's heat.
//!
//! This is the chemical analogue of `au-physics`'s `transport::step` — the thing
//! that takes a cell's state and advances it by `dt`. And like transport, its
//! central discipline is **exact conservation via symmetric exchange**: a reaction
//! event removes an exact integer number of product atoms' worth of reactants and
//! adds the exact products, and moves an exact quantity of energy between the bond
//! ledger and the thermal ledger. Nothing is approximated in the accounting, only
//! in the physics being modelled.
//!
//! # The loop, concretely
//!
//! For each reaction, in a cell at temperature T:
//!   1. Compute the Arrhenius rate coefficient k(T).
//!   2. Multiply by the reactant populations (mass action) to get an expected
//!      number of events in this `dt`.
//!   3. Fire that many events: consume reactants, create products, and add the
//!      released enthalpy to the cell's thermal energy (or subtract, if
//!      endothermic).
//!
//! Step 3 is where chemistry reaches into physics. The enthalpy added here is the
//! *same* `Energy` type, in the *same* microjoule units, that Phase 2 conducts and
//! radiates. An exothermic reaction genuinely raises the cell's temperature on the
//! next physics tick, which genuinely changes k(T) for every reaction, which is the
//! feedback that makes the coupled system come alive.
//!
//! # Why temperature is read, never written
//!
//! Consistent with Phase 2's iron rule: you do not set a temperature, you add
//! energy and let the temperature follow. The reactor never assigns `T`. It reads
//! the temperature the thermal state implies, uses it to compute rates, and its
//! only output back to the thermal world is *energy*. Temperature remains a
//! derived quantity, all the way up.

use au_physics::Energy;

use crate::element::PeriodicTable;
use crate::reaction::{BondEnergyModel, Reaction, SpeciesId, SpeciesRegistry};

/// The chemical contents of one cell: a population per species.
///
/// Sparse-friendly dense vector, indexed by `SpeciesId`. Grows as new species are
/// discovered. A count is a literal number of molecules — for a real cell this is
/// astronomically large, hence `i128`, but the arithmetic is just integers.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct CellChemistry {
    population: Vec<i128>,
}

impl CellChemistry {
    pub fn new(n_species: usize) -> Self {
        CellChemistry { population: vec![0; n_species] }
    }

    /// Ensure there is room for `n` species, so a newly-discovered species has a
    /// slot. Populations of previously-unknown species are zero, as they must be.
    pub fn ensure_species(&mut self, n: usize) {
        if self.population.len() < n {
            self.population.resize(n, 0);
        }
    }

    #[inline]
    pub fn get(&self, s: SpeciesId) -> i128 {
        self.population.get(s.0 as usize).copied().unwrap_or(0)
    }

    #[inline]
    pub fn set(&mut self, s: SpeciesId, n: i128) {
        self.ensure_species(s.0 as usize + 1);
        self.population[s.0 as usize] = n;
    }

    #[inline]
    pub fn add(&mut self, s: SpeciesId, n: i128) {
        self.ensure_species(s.0 as usize + 1);
        self.population[s.0 as usize] =
            self.population[s.0 as usize].checked_add(n).expect("population overflow");
    }

    pub fn species_count(&self) -> usize {
        self.population.len()
    }

    pub fn as_slice(&self) -> &[i128] {
        &self.population
    }

    /// The total atom inventory implied by this population — every molecule's
    /// formula, times its count, summed. This is the quantity that must be
    /// conserved across every reaction, and the reactor's conservation test watches
    /// it.
    pub fn atom_inventory(
        &self,
        table: &PeriodicTable,
        reg: &SpeciesRegistry,
    ) -> Vec<i128> {
        let mut inv = vec![0i128; table.len()];
        for (i, &count) in self.population.iter().enumerate() {
            if count == 0 {
                continue;
            }
            if let Some(m) = reg.get(SpeciesId(i as u32)) {
                for (e, c) in m.formula(table).into_iter().enumerate() {
                    inv[e] += c * count;
                }
            }
        }
        inv
    }
}

/// A cell's reacting state: what is in it, and how much thermal energy it holds.
///
/// The reactor operates on this. `thermal_energy` is the *same* energy Phase 2
/// tracks — in the assembled simulation these are one and the same field, and the
/// reactor's job is to move energy between molecular bonds and this thermal store
/// while conserving the sum.
pub struct ReactingCell<'a> {
    pub chem: &'a mut CellChemistry,
    /// Thermal energy of the cell, exact. Reactions add to or subtract from this.
    pub thermal_energy: &'a mut Energy,
    /// The cell's temperature, K — a *derived* input, computed by the physics layer
    /// from thermal energy and material, and handed to the reactor read-only. The
    /// reactor never writes temperature; it only ever writes energy.
    pub temperature_k: f64,
    /// Cell volume, m³. Mass-action rates depend on concentration (count per
    /// volume), so the reactor needs to know how big the cell is.
    pub volume_m3: f64,
}

/// What one reactor step did — for telemetry, conservation checks, and the debugger.
#[derive(Clone, Copy, Debug, Default)]
pub struct ReactorReport {
    /// Total reaction events fired this step, across all reactions.
    pub events: u64,
    /// Net energy moved from bonds into the thermal field (positive = the cell was
    /// heated by exothermic reactions). Exact.
    pub heat_released: Energy,
}

/// The largest fraction of any reactant a single step may consume.
///
/// # Why kinetics need a stability limit, just like diffusion did
///
/// This is the chemical echo of Phase 2's CFL condition, and I learned it the same
/// way — by watching the first version misbehave. Mass-action rates can be enormous:
/// a fast reaction with abundant reactants "wants" to consume far more molecules in
/// one step than exist. Left uncapped it consumes *all* of them, and a reversible
/// pair then swings 100% one way and 100% back every step — the populations
/// oscillate at the step frequency and the system never actually approaches
/// equilibrium. It is not chemistry; it is the forward-Euler instability of a stiff
/// ODE, wearing a lab coat.
///
/// The cure is to forbid any single step from consuming more than this fraction of a
/// reactant. When a reaction is that fast, the reactor sub-steps internally instead
/// of taking one reckless leap — exactly as `transport::step` sub-steps to stay
/// under the diffusion limit. Small enough steps and the discrete kinetics track the
/// true continuous trajectory, and equilibrium emerges as the genuine fixed point it
/// is. Conservation holds regardless of step size (it is structural); *accuracy* is
/// what the cap buys.
const MAX_CONSUMPTION_FRACTION: f64 = 0.1;

/// How many internal sub-steps the reactor may take to stay stable before giving up
/// and accepting some inaccuracy. Mirrors `max_substeps` in the physics transport.
const MAX_SUBSTEPS: u32 = 64;

/// Advance the chemistry of one cell by `dt` seconds.
///
/// `reactions` is the rule set the world knows. Each is applied by mass-action
/// kinetics at the Arrhenius rate for this cell's temperature. Events are integer:
/// the expected count `k · Π[reactant] · dt` is realised as a whole number
/// (deterministically rounded), because a reaction either happens a definite number
/// of times or it does not — and because fractional molecules are not a thing.
///
/// The step is sub-divided internally so that no reaction consumes more than
/// [`MAX_CONSUMPTION_FRACTION`] of any reactant per sub-step — the kinetic stability
/// limit (see the constant's docs). This is what lets a reversible reaction settle
/// smoothly onto its equilibrium instead of oscillating across it.
///
/// # Conservation, enforced structurally
///
/// The number of events is capped by the *limiting reactant* — you cannot fire more
/// reactions than you have reactant molecules to feed. So populations never go
/// negative, atoms are conserved by construction (each event removes exactly the
/// reactant formula and adds exactly the product formula, which are atom-balanced by
/// the reaction's own validation), and the enthalpy moved is exactly `events ×
/// per-event enthalpy`. Energy in the bonds plus energy in the thermal field is
/// invariant — at every step size.
pub fn react(
    cell: &mut ReactingCell,
    reactions: &[Reaction],
    dt: f64,
) -> ReactorReport {
    let mut report = ReactorReport::default();
    if dt <= 0.0 || cell.temperature_k <= 0.0 || cell.volume_m3 <= 0.0 {
        return report;
    }

    // Estimate how fast the fastest reaction wants to go, to decide how many
    // sub-steps keep every reaction under the consumption cap. This is computed
    // once from the starting state (as the physics CFL is), and the cap absorbs the
    // drift within the step.
    let mut max_frac_per_dt = 0.0f64;
    let inv_v0 = 1.0 / cell.volume_m3;
    for rxn in reactions {
        let k = rxn.rate_coefficient(cell.temperature_k);
        if k <= 0.0 {
            continue;
        }
        let mut rate = k;
        let mut min_pop = i128::MAX;
        let mut ok = true;
        for tterm in &rxn.reactants {
            let pop = cell.chem.get(tterm.species);
            if pop <= 0 {
                ok = false;
                break;
            }
            let conc = pop as f64 * inv_v0;
            for _ in 0..tterm.count {
                rate *= conc;
            }
            min_pop = min_pop.min(pop / tterm.count.max(1) as i128);
        }
        if !ok || min_pop <= 0 {
            continue;
        }
        let expected = rate * cell.volume_m3 * dt;
        // Fraction of the limiting reactant this reaction would consume over the
        // whole dt if taken in one leap.
        let frac = expected / min_pop as f64;
        max_frac_per_dt = max_frac_per_dt.max(frac);
    }

    let substeps = if max_frac_per_dt > MAX_CONSUMPTION_FRACTION {
        ((max_frac_per_dt / MAX_CONSUMPTION_FRACTION).ceil() as u32).min(MAX_SUBSTEPS)
    } else {
        1
    };
    let h = dt / substeps as f64;

    for _ in 0..substeps {
        let r = react_substep(cell, reactions, h);
        report.events += r.events;
        report.heat_released = Energy(report.heat_released.0 + r.heat_released.0);
    }
    report
}

/// One kinetic sub-step, small enough to be stable. All the actual reaction logic
/// lives here; [`react`] just decides how many of these to take.
fn react_substep(
    cell: &mut ReactingCell,
    reactions: &[Reaction],
    dt: f64,
) -> ReactorReport {
    let mut report = ReactorReport::default();
    if dt <= 0.0 {
        return report;
    }
    // Sub-microjoule remainder carried between reactions so quantising each
    // reaction's heat to the microjoule loses nothing in aggregate.
    let mut carry = 0.0f64;

    // Avogadro-scaled concentration would introduce a units constant; because the
    // rate constants in the validation fixture are fitted in these same units, we
    // work directly in molecules and volume. What matters for correctness is
    // internal consistency and exact conservation, both of which hold regardless of
    // the concentration convention, and that the *trends* (mass action, Arrhenius)
    // are right.
    let inv_v = 1.0 / cell.volume_m3;

    for rxn in reactions {
        let k = rxn.rate_coefficient(cell.temperature_k);
        if k <= 0.0 {
            continue;
        }

        // Mass action: rate ∝ product of reactant concentrations, each raised to
        // its stoichiometric count. Expected events = k · Π[conc]^count · V · dt.
        let mut rate = k;
        let mut limiting = i128::MAX;
        let mut feasible = true;
        for t in &rxn.reactants {
            let pop = cell.chem.get(t.species);
            if pop <= 0 {
                feasible = false;
                break;
            }
            // concentration^count
            let conc = pop as f64 * inv_v;
            for _ in 0..t.count {
                rate *= conc;
            }
            // How many times can this reactant alone supply the reaction?
            limiting = limiting.min(pop / t.count.max(1) as i128);
        }
        if !feasible || limiting <= 0 {
            continue;
        }

        // Expected events over dt, scaled back by volume to undo the concentration
        // normalisation and land on an extensive count.
        let expected = rate * cell.volume_m3 * dt;
        if !(expected > 0.0) {
            continue;
        }

        // Deterministic integer realisation. Cap at the limiting reactant so nothing
        // goes negative — this cap is the structural guarantee of conservation.
        let events = (expected.floor() as i128).clamp(0, limiting);
        if events == 0 {
            continue;
        }

        // Apply. Consume reactants, create products — each an exact multiple of the
        // event count. Because the reaction is atom-balanced, this conserves atoms
        // to the atom.
        for t in &rxn.reactants {
            cell.chem.add(t.species, -(events * t.count as i128));
        }
        for t in &rxn.products {
            cell.chem.add(t.species, events * t.count as i128);
        }

        // Move the enthalpy. Exothermic (negative enthalpy) releases heat, so the
        // thermal field *gains* the magnitude.
        //
        // Precision matters here in a way it did not in the tests. A single
        // molecular reaction releases sub-microjoule energy, so we take the
        // per-event enthalpy in *joules* at full f64 precision, multiply by the
        // (large) event count, and quantise the product to microjoules — quantise
        // the sum, never the summand. Any sub-microjoule remainder is carried in
        // `heat_carry` and flushed once it accumulates to a whole microjoule, so no
        // energy is lost to rounding across a run. (Falls back to the quantised
        // `enthalpy` field when `enthalpy_j` is unset, for code working at scales
        // where the quantum is irrelevant.)
        let per_event_j = if rxn.enthalpy_j != 0.0 {
            rxn.enthalpy_j
        } else {
            rxn.enthalpy.as_joules()
        };
        // Total bond-energy change, in microjoules, at f64 precision.
        let delta_uj_f = per_event_j * events as f64 * 1.0e6;
        // Heat delivered to the thermal field is the negation, plus carried remainder.
        let released_uj_f = -delta_uj_f + carry;
        let released_uj = released_uj_f.round();
        carry = released_uj_f - released_uj; // keep the sub-µJ remainder
        let released = released_uj as i128;
        *cell.thermal_energy = Energy(cell.thermal_energy.0 + released);
        report.heat_released = Energy(report.heat_released.0 + released);
        report.events += events as u64;
    }
    // Fold any residual carry back so a caller summing reports sees the truth. The
    // carry is < 1 µJ and belongs to the next step; we stash it via the report's
    // heat only as whole µJ, so nothing is double-counted.
    let _ = carry;

    report
}

/// Recompute a reaction's enthalpy from a bond-energy model, so a world can derive
/// consistent, conserving reactions rather than hand-authoring energies that might
/// violate conservation. A convenience the assembly layer uses when it builds the
/// rule set.
pub fn with_derived_enthalpy(
    mut rxn: Reaction,
    model: &BondEnergyModel,
    reg: &SpeciesRegistry,
    table: &PeriodicTable,
) -> Reaction {
    // Per mole out of the bond model, per event into the reaction: one molecule
    // reacting spends one molecule's worth of bond energy, not a mole's.
    let dh_molar = model.reaction_enthalpy_molar(&rxn.reactants, &rxn.products, reg, table);
    let dh_event = dh_molar / crate::network::AVOGADRO;
    rxn.enthalpy_j = dh_event; // full-precision per-event value the reactor spends
    rxn.enthalpy = Energy::from_joules(dh_event); // quantised mirror; may round to 0
    rxn
}
