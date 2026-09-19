//! The world's declared chemistry — its elements, its molecules, its reactions —
//! built from config at boot, deterministically.
//!
//! This is the chemical analogue of `MaterialRegistry::from_kv`, and it exists for
//! the same reason: a reloaded world must run on *exactly* the chemistry it was
//! computed with. Because the whole vocabulary is a pure function of config, it
//! need not be serialized into snapshots at all — it is rebuilt bit-identically on
//! load, the way the material table is. What *is* serialized is the per-cell
//! populations (through the chemistry columns); the meaning of those populations —
//! which id is which molecule — is reconstructed from config.
//!
//! # What "declaring chemistry" means, and does not mean
//!
//! Config names the elements that exist, the molecules that *can* exist, and the
//! reaction rules that *may* fire. That is the vocabulary. It is emphatically not
//! the outcome: nothing here says a reaction *will* fire, or how far, or how hot a
//! cell gets, or which molecule dominates. Those are decided by the reactor from
//! populations, temperature, and rate laws — emergent, every time. Declaring that
//! carbon and water may exist is the same kind of act as declaring which elements
//! the periodic table contains. It is the opposite of scripting a wolf.
//!
//! # The config schema
//!
//! ```text
//! chem.enabled = true
//!
//! chem.element.count = 2
//! chem.element.0.z = 1
//! chem.element.0.symbol = H
//! chem.element.0.mass_mda = 1008
//! chem.element.0.valence = 1
//! chem.element.0.electronegativity_c = 220
//! chem.element.1.z = 8   (… O …)
//!
//! chem.species.count = 3
//! chem.species.0.atoms = 1,1          # two hydrogen atoms (by Z)
//! chem.species.0.bonds = 0-1:1        # atom0–atom1, single bond   → H2
//! chem.species.1.atoms = 8,8
//! chem.species.1.bonds = 0-1:1        # O2
//! chem.species.2.atoms = 1,8
//! chem.species.2.bonds = 0-1:1        # H–O
//!
//! chem.reaction.count = 2
//! chem.reaction.0.reactants = 0:1,1:1     # 1×species0 + 1×species1
//! chem.reaction.0.products  = 2:2         # → 2×species2
//! chem.reaction.0.activation_j = 40000
//! chem.reaction.0.pre_exponential = 1000
//! chem.reaction.0.enthalpy_j = -4.98e-20  # per-event; if absent, derived from bonds
//! chem.reaction.1. … (the reverse) …
//!
//! chem.bond.base_per_order = 200000       # optional; defaults match au-chem
//! chem.bond.electroneg_coefficient = 100000
//! ```
//!
//! Bond edges are `atomA-atomB:order`, order ∈ {1,2,3}. Reaction terms are
//! `speciesIndex:count`, comma-separated.

use std::collections::BTreeMap;

use au_chem::{
    Atom, Bond, BondEnergyModel, BondOrder, Element, Molecule, PeriodicTable, Reaction,
    SpeciesId, SpeciesRegistry, Term, Z,
};
use au_physics::Energy;

/// A world's complete chemical vocabulary. Immutable after boot.
pub struct Chemistry {
    pub table: PeriodicTable,
    pub species: SpeciesRegistry,
    pub reactions: Vec<Reaction>,
    pub bond_model: BondEnergyModel,
    /// Cell volume, m³, for mass-action concentrations. One value for the whole
    /// world in Phase 3b (a per-cell volume is a Phase 4 refinement, when the
    /// planet gives cells real sizes).
    pub cell_volume_m3: f64,

    /// Open-ended mode: reactions are *generated* from molecular structure and
    /// new molecules are discovered at runtime (Phase 5a). `None` means the
    /// declared Phase 3 chemistry, unchanged. The declared species become the
    /// **seeds** — usually bare atoms — and everything else must be reached.
    pub open: Option<au_chem::NetworkRules>,
    /// Diffusion coefficient, m²/s, of a **1-dalton reference particle** in the
    /// host fluid; 0 disables species transport. Every species then derives its
    /// own rate by Graham's scaling, D_s = D₀·√(1 Da / m_s) — lighter molecules
    /// move faster — with the mass read off the discovered molecular graph. The
    /// *shape* of transport is derived from structure; only the overall scale
    /// is declared, exactly as with the reaction barrier.
    pub diffusion_m2_s: f64,

    /// Derived catalysis (Phase 5c). `None` — the default — means no catalysis
    /// at all, and a world configured that way runs byte-identically to Phase
    /// 5b. When enabled, molecules that can grip two reactants at once lower
    /// that reaction's barrier in *both* directions, and the catalysed variants
    /// join the reaction set with the catalyst written on both sides. That is
    /// what makes `A + B + AB → 2 AB` expressible: the first stoichiometry in
    /// this engine's life in which a species can increase its own number.
    pub catalysis: Option<au_chem::CatalysisRules>,

    /// Ceiling on how many catalysed variants a world will carry. Catalysis is
    /// quadratic — every reaction × every candidate catalyst — and an unbounded
    /// set would make the reactor's per-cell loop unaffordable. When the
    /// derivation overflows this, the *strongest* catalysts are kept (largest
    /// barrier reduction first, ties broken by id) and the rest dropped. A
    /// declared edge of the world, like the species pool: stated, deterministic,
    /// and reported rather than hit silently.
    pub max_catalysed: usize,

    /// Self-assembly. `None` — the default — means no membranes, and a world
    /// configured that way runs byte-identically to Phase 5e. When enabled,
    /// amphiphilic species past their critical concentration become surface,
    /// and that surface slows every exchange across the cell's faces. Nothing is
    /// stored: the membrane is a function of the populations, recomputed each
    /// tick exactly as temperature is recomputed from energy.
    pub membrane: Option<au_chem::MembraneRules>,

    /// Open boundaries for matter (Phase 5h). Empty — the default — means a
    /// sealed world, and a world configured that way runs byte-identically to
    /// Phase 5g. Each entry is a cell index within every active chunk and the
    /// port that holds it: a source at a fixed composition, or a sink at zero.
    ///
    /// Ports live in config rather than in the snapshot for the same reason
    /// every other rule does: they are a statement about the *world's* setup,
    /// not about its state. A resumed world gets them from the config it is
    /// resumed with, exactly as it gets its chemistry and its materials. What
    /// the snapshot must carry is the *ledger* — what has already crossed —
    /// because that is history.
    pub ports: Vec<(usize, au_chem::Port)>,

    /// Column capacity — the world's shelf space for species. In closed mode it
    /// equals the declared species count; in open mode it is `chem.max_species`,
    /// and discovery past it is recorded but never stocked (see the system).
    pub pool: usize,
}

impl Chemistry {
    pub fn species_count(&self) -> usize {
        self.species.len()
    }

    /// How many population columns this world carries per cell.
    pub fn column_count(&self) -> usize {
        self.pool
    }

    /// Build the vocabulary from config. Returns `None` if chemistry is disabled,
    /// `Err` if the config is malformed — a malformed chemistry is a hard error,
    /// not a shrug, because every cell's populations would otherwise be meaningless.
    pub fn from_kv(kv: &BTreeMap<String, String>) -> Result<Option<Chemistry>, String> {
        let enabled = kv.get("chem.enabled").map(|s| s == "true").unwrap_or(false);
        if !enabled {
            return Ok(None);
        }

        // ── Elements. ───────────────────────────────────────────────────────
        let mut table = PeriodicTable::new();
        let n_elem: usize = parse_at(kv, "chem.element.count").unwrap_or(Ok(0))?;
        for i in 0..n_elem {
            let z = Z(parse_at(kv, &format!("chem.element.{}.z", i))
                .ok_or_else(|| format!("missing chem.element.{}.z", i))??);
            let symbol = kv
                .get(&format!("chem.element.{}.symbol", i))
                .cloned()
                .unwrap_or_else(|| format!("E{}", z.0));
            // Symbols must outlive the run; leak the few element symbols so the
            // `&'static str` contract of `Element` is satisfied. A handful of
            // short strings, once, at boot.
            let symbol: &'static str = Box::leak(symbol.into_boxed_str());
            table.add(Element {
                z,
                symbol,
                mass_mda: parse_at(kv, &format!("chem.element.{}.mass_mda", i))
                    .ok_or_else(|| format!("missing chem.element.{}.mass_mda", i))??,
                valence: parse_at(kv, &format!("chem.element.{}.valence", i))
                    .ok_or_else(|| format!("missing chem.element.{}.valence", i))??,
                electronegativity_c: parse_at(
                    kv,
                    &format!("chem.element.{}.electronegativity_c", i),
                )
                .ok_or_else(|| format!("missing chem.element.{}.electronegativity_c", i))??,
            });
        }

        // ── Bond model (optional; defaults match au-chem's). ────────────────
        let bond_model = BondEnergyModel {
            base_per_order: kv
                .get("chem.bond.base_per_order")
                .and_then(|s| s.parse().ok())
                .unwrap_or(BondEnergyModel::default().base_per_order),
            electroneg_coefficient: kv
                .get("chem.bond.electroneg_coefficient")
                .and_then(|s| s.parse().ok())
                .unwrap_or(BondEnergyModel::default().electroneg_coefficient),
        };

        // ── Species (molecules as graphs). Interned in declared order, so the
        //    ids are stable across runs and across save/resume. ───────────────
        let mut species = SpeciesRegistry::new();
        let n_species: usize = parse_at(kv, "chem.species.count").unwrap_or(Ok(0))?;
        for i in 0..n_species {
            let mol = parse_molecule(kv, i)?;
            let id = species.intern(mol);
            // Interning in order must yield sequential ids; if a duplicate molecule
            // is declared it would collapse two indices and corrupt the column
            // mapping. Refuse.
            if id.0 as usize != i {
                return Err(format!(
                    "chem.species.{} is a duplicate of species {} — every declared species must be distinct",
                    i, id.0
                ));
            }
        }

        // ── Reactions. ──────────────────────────────────────────────────────
        let mut reactions = Vec::new();
        let n_rxn: usize = parse_at(kv, "chem.reaction.count").unwrap_or(Ok(0))?;
        for i in 0..n_rxn {
            reactions.push(parse_reaction(kv, i, &species, &table, &bond_model)?);
        }

        // Every reaction must balance atoms. A reaction that does not is a hole in
        // conservation, and the whole tower rests on there being none.
        for (i, r) in reactions.iter().enumerate() {
            if !r.is_atom_balanced(&table, &species) {
                return Err(format!(
                    "chem.reaction.{} does not balance atoms — it creates or destroys matter",
                    i
                ));
            }
        }

        let cell_volume_m3 = kv
            .get("chem.cell_volume_m3")
            .and_then(|s| s.parse().ok())
            .unwrap_or(1.0);

        // ── Open-ended chemistry (Phase 5a). The declared species above become
        //    seeds; the network generator does the rest. Declared *reactions*
        //    still work alongside generated ones is a combination nobody needs
        //    yet, so it is refused loudly rather than half-supported.
        let open_ended = kv.get("chem.open_ended").map(|s| s == "true").unwrap_or(false);
        let (open, pool) = if open_ended {
            if !reactions.is_empty() {
                return Err(
                    "chem.open_ended worlds derive their reactions; declaring chem.reaction.* \
                     alongside is refused rather than half-supported"
                        .into(),
                );
            }
            // The old dial is gone, and it must fail loudly rather than be
            // silently reinterpreted: `chem.pre_exponential` was a single fitted
            // number applied to unimolecular and bimolecular reactions alike,
            // which cannot be dimensionally correct for both. Prefactors are now
            // derived per reaction from molecular mass and size, and the only
            // thing left to declare is the dimensionless steric factor.
            if kv.contains_key("chem.pre_exponential") {
                return Err("chem.pre_exponential no longer exists: reaction prefactors are \
                            derived from molecular structure (collision theory). The only \
                            declared kinetic scale is chem.steric_factor (dimensionless, \
                            default 1.0). See au_chem::kinetics."
                    .to_string());
            }
            let rules = au_chem::NetworkRules {
                intrinsic_barrier_j_mol: kv
                    .get("chem.barrier_j_mol")
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(60_000.0),
                kinetics: au_chem::KineticRules {
                    steric_factor: kv
                        .get("chem.steric_factor")
                        .and_then(|s| s.parse().ok())
                        .unwrap_or(1.0),
                    atom_radius_m: kv
                        .get("chem.atom_radius_m")
                        .and_then(|s| s.parse().ok())
                        .unwrap_or(1.8e-10),
                    encounter_volume_m3: kv
                        .get("chem.encounter_volume_m3")
                        .and_then(|s| s.parse().ok())
                        .unwrap_or(1.0e-28),
                },
                max_atoms: kv
                    .get("chem.max_atoms")
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(8),
            };
            let pool: usize = kv
                .get("chem.max_species")
                .and_then(|s| s.parse().ok())
                .unwrap_or(64);
            if pool < species.len() {
                return Err(format!(
                    "chem.max_species = {} is smaller than the {} declared seed species",
                    pool,
                    species.len()
                ));
            }
            (Some(rules), pool)
        } else {
            (None, species.len())
        };

        let diffusion_m2_s = kv
            .get("chem.diffusion_m2_s")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        // Catalysis: off unless asked for, so every earlier world is untouched.
        let catalysis = if kv.get("chem.catalysis").map(|v| v == "true").unwrap_or(false) {
            let d = au_chem::CatalysisRules::default();
            Some(au_chem::CatalysisRules {
                max_reduction_frac: kv
                    .get("chem.catalysis.max_reduction_frac")
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(d.max_reduction_frac),
                min_grip_j_mol: kv
                    .get("chem.catalysis.min_grip_j_mol")
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(d.min_grip_j_mol),
                grip_half_j_mol: kv
                    .get("chem.catalysis.grip_half_j_mol")
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(d.grip_half_j_mol),
            })
        } else {
            None
        };
        let membrane = if kv.get("chem.membrane").map(|v| v == "true").unwrap_or(false) {
            let d = au_chem::MembraneRules::default();
            Some(au_chem::MembraneRules {
                cmc_ref_per_m3: kv.get("chem.membrane.cmc_ref_per_m3")
                    .and_then(|s| s.parse().ok()).unwrap_or(d.cmc_ref_per_m3),
                cmc_decade: kv.get("chem.membrane.cmc_decade")
                    .and_then(|s| s.parse().ok()).unwrap_or(d.cmc_decade),
                min_score: kv.get("chem.membrane.min_score")
                    .and_then(|s| s.parse().ok()).unwrap_or(d.min_score),
                area_per_molecule_m2: kv.get("chem.membrane.area_per_molecule_m2")
                    .and_then(|s| s.parse().ok()).unwrap_or(d.area_per_molecule_m2),
                residual_permeability: kv.get("chem.membrane.residual_permeability")
                    .and_then(|s| s.parse().ok()).unwrap_or(d.residual_permeability),
            })
        } else {
            None
        };

        let max_catalysed = kv
            .get("chem.catalysis.max_variants")
            .and_then(|s| s.parse().ok())
            .unwrap_or(4096);
        if catalysis.is_some() && open.is_none() {
            return Err(
                "chem.catalysis requires chem.open_ended: catalysis is derived from molecular \
                 structure, and a declared reaction list has no structure to derive it from"
                    .to_string(),
            );
        }

        let ports = parse_ports(kv, &species)?;

        Ok(Some(Chemistry {
            table,
            species,
            reactions,
            bond_model,
            cell_volume_m3,
            open,
            pool,
            diffusion_m2_s,
            catalysis,
            max_catalysed,
            membrane,
            ports,
        }))
    }
}

/// Parse the port declarations.
///
/// ```text
///   chem.reservoir.source = "0:1@200000,2@200000; 5:1@50"
///   chem.reservoir.sink   = "40; 41"
/// ```
///
/// Cells are indices within a chunk; species are ids. Ids rather than formulae
/// because in open mode most species do not exist yet at boot — the world has
/// not discovered them — and the ones a port can sensibly name are the declared
/// seeds, whose ids are their declaration order.
///
/// A malformed port is a hard error rather than a shrug. Silently ignoring an
/// unparseable boundary condition would give a world that *looks* open, runs to
/// equilibrium, and tells nobody why.
fn parse_ports(
    kv: &BTreeMap<String, String>,
    species: &au_chem::SpeciesRegistry,
) -> Result<Vec<(usize, au_chem::Port)>, String> {
    let mut out: Vec<(usize, au_chem::Port)> = Vec::new();

    if let Some(spec) = kv.get("chem.reservoir.source") {
        for group in spec.split(';').map(str::trim).filter(|g| !g.is_empty()) {
            let (cell_s, rest) = group
                .split_once(':')
                .ok_or_else(|| format!("chem.reservoir.source: expected `cell:sp@pop` in `{}`", group))?;
            let cell: usize = cell_s
                .trim()
                .parse()
                .map_err(|_| format!("chem.reservoir.source: bad cell index `{}`", cell_s))?;
            let mut hold = Vec::new();
            for item in rest.split(',').map(str::trim).filter(|i| !i.is_empty()) {
                let (sp, pop) = item
                    .split_once('@')
                    .ok_or_else(|| format!("chem.reservoir.source: expected `sp@pop` in `{}`", item))?;
                let sp: u32 = sp
                    .trim()
                    .parse()
                    .map_err(|_| format!("chem.reservoir.source: bad species id `{}`", sp))?;
                let pop: i128 = pop
                    .trim()
                    .parse()
                    .map_err(|_| format!("chem.reservoir.source: bad population `{}`", pop))?;
                if pop < 0 {
                    return Err("chem.reservoir.source: a held population cannot be negative".into());
                }
                if species.get(au_chem::SpeciesId(sp)).is_none() {
                    return Err(format!(
                        "chem.reservoir.source: species {} is not declared; a port may only \
                         hold a species the world knows about at boot",
                        sp
                    ));
                }
                hold.push((au_chem::SpeciesId(sp), pop));
            }
            hold.sort_by_key(|(s, _)| s.0);
            out.push((cell, au_chem::Port::Source { hold }));
        }
    }

    if let Some(spec) = kv.get("chem.reservoir.sink") {
        for item in spec.split(';').map(str::trim).filter(|i| !i.is_empty()) {
            let cell: usize = item
                .parse()
                .map_err(|_| format!("chem.reservoir.sink: bad cell index `{}`", item))?;
            out.push((cell, au_chem::Port::Sink));
        }
    }

    // A cell cannot be both a source and a sink: the two would fight, and which
    // one won would depend on declaration order — a rule whose result depends on
    // the order it was written in is not a rule.
    out.sort_by_key(|(c, _)| *c);
    for w in out.windows(2) {
        if w[0].0 == w[1].0 {
            return Err(format!(
                "chem.reservoir: cell {} is declared as a port twice; a cell has one boundary \
                 condition or none",
                w[0].0
            ));
        }
    }
    Ok(out)
}

fn parse_at<T: std::str::FromStr>(
    kv: &BTreeMap<String, String>,
    key: &str,
) -> Option<Result<T, String>> {
    kv.get(key).map(|s| s.parse::<T>().map_err(|_| format!("{} is not valid", key)))
}

/// Parse `chem.species.{i}` into a molecule graph.
fn parse_molecule(kv: &BTreeMap<String, String>, i: usize) -> Result<Molecule, String> {
    let atoms_s = kv
        .get(&format!("chem.species.{}.atoms", i))
        .ok_or_else(|| format!("missing chem.species.{}.atoms", i))?;
    let atoms: Vec<Atom> = atoms_s
        .split(',')
        .map(|t| {
            t.trim()
                .parse::<u8>()
                .map(|z| Atom { z: Z(z) })
                .map_err(|_| format!("chem.species.{}.atoms has a bad atomic number '{}'", i, t))
        })
        .collect::<Result<_, _>>()?;

    let bonds_s = kv.get(&format!("chem.species.{}.bonds", i)).cloned().unwrap_or_default();
    let mut bonds = Vec::new();
    if !bonds_s.trim().is_empty() {
        for edge in bonds_s.split(',') {
            let edge = edge.trim();
            // Format "a-b:order".
            let (pair, order) = edge
                .split_once(':')
                .ok_or_else(|| format!("chem.species.{}.bonds edge '{}' missing :order", i, edge))?;
            let (a, b) = pair
                .split_once('-')
                .ok_or_else(|| format!("chem.species.{}.bonds edge '{}' missing a-b", i, edge))?;
            let a: u16 = a.trim().parse().map_err(|_| format!("bad atom index '{}'", a))?;
            let b: u16 = b.trim().parse().map_err(|_| format!("bad atom index '{}'", b))?;
            let o: u8 = order.trim().parse().map_err(|_| format!("bad bond order '{}'", order))?;
            let order = BondOrder::from_u8(o)
                .ok_or_else(|| format!("chem.species.{}.bonds order {} not in 1..3", i, o))?;
            if a as usize >= atoms.len() || b as usize >= atoms.len() {
                return Err(format!(
                    "chem.species.{}.bonds edge '{}' references an atom outside the molecule",
                    i, edge
                ));
            }
            bonds.push(Bond::new(a, b, order));
        }
    }

    Ok(Molecule::new(atoms, bonds))
}

/// Parse `chem.reaction.{i}` into a reaction rule.
fn parse_reaction(
    kv: &BTreeMap<String, String>,
    i: usize,
    species: &SpeciesRegistry,
    table: &PeriodicTable,
    bond_model: &BondEnergyModel,
) -> Result<Reaction, String> {
    let terms = |field: &str| -> Result<Vec<Term>, String> {
        let s = kv
            .get(&format!("chem.reaction.{}.{}", i, field))
            .ok_or_else(|| format!("missing chem.reaction.{}.{}", i, field))?;
        s.split(',')
            .map(|t| {
                let t = t.trim();
                let (idx, count) = t
                    .split_once(':')
                    .ok_or_else(|| format!("reaction {} term '{}' must be species:count", i, t))?;
                let idx: u32 = idx.trim().parse().map_err(|_| format!("bad species index '{}'", idx))?;
                let count: u32 = count.trim().parse().map_err(|_| format!("bad count '{}'", count))?;
                if idx as usize >= species.len() {
                    return Err(format!("reaction {} references species {} which is not declared", i, idx));
                }
                Ok(Term { species: SpeciesId(idx), count })
            })
            .collect()
    };

    let reactants = terms("reactants")?;
    let products = terms("products")?;

    let activation_j = parse_at(kv, &format!("chem.reaction.{}.activation_j", i))
        .ok_or_else(|| format!("missing chem.reaction.{}.activation_j", i))??;
    let pre_exponential = parse_at(kv, &format!("chem.reaction.{}.pre_exponential", i))
        .ok_or_else(|| format!("missing chem.reaction.{}.pre_exponential", i))??;

    // Enthalpy: either given explicitly (per event, joules) or derived from bond
    // energies. Deriving is the safer default — a derived enthalpy cannot violate
    // energy conservation, because it is *defined* as the bond-energy difference.
    let enthalpy_j: f64 = match kv.get(&format!("chem.reaction.{}.enthalpy_j", i)) {
        Some(s) => s.parse().map_err(|_| format!("bad enthalpy_j for reaction {}", i))?,
        None => {
            // Bond energies are tabulated per mole; a reaction event consumes one
            // molecule's share. Omitting this division is what drove a 1e-15 kg
            // cell to absolute zero in Phase 5i step 2.
            bond_model.reaction_enthalpy_molar(&reactants, &products, species, table)
                / au_chem::AVOGADRO
        }
    };

    Ok(Reaction {
        reactants,
        products,
        activation_j,
        enthalpy: Energy::from_joules(enthalpy_j),
        enthalpy_j,
        pre_exponential,
    })
}
