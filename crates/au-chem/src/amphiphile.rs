//! Amphiphilicity — how a molecule comes to have two natures at once.
//!
//! # Why a network needs this
//!
//! Phase 5c found sets of reactions that collectively make themselves. What such
//! a set does not have is an *inside*. A RAF is a property of the whole pot: one
//! soup, so nothing can be an individual, and with no individuals there is
//! nothing for selection to act on. Every organism that has ever existed solved
//! this the same way — with a bag.
//!
//! Bags are not a new kind of physics. They follow from a molecule having a
//! polar end and an apolar end: water can order itself around the polar parts
//! and cannot around the others, so the cheapest arrangement hides the apolar
//! parts from the water. Amphiphiles therefore assemble into films, micelles and
//! bilayers *by themselves*, with no machinery — which is exactly the sort of
//! thing this project is supposed to let happen rather than script.
//!
//! So the question here is narrow and structural: **given a molecular graph the
//! world discovered, can it be cut into a polar piece and an apolar piece?**
//! Nothing in this module builds a membrane. It measures a property; a later
//! layer decides what a cell full of such molecules does.
//!
//! # The measure
//!
//! Every bond between unlike elements is polar in proportion to the difference
//! in electronegativity — the same quantity, from the same table, that prices
//! every bond in [`crate::reaction`]. An atom's polarity is the sum over its
//! bonds of `order × |Δ electronegativity|`, so an oxygen holding a hydrogen is
//! polar and a carbon surrounded by carbons and hydrogens is not.
//!
//! The molecule is then scanned at every **bridge** bond — every bond whose
//! removal separates it in two. Rings have none, and so correctly score zero: a
//! ring has no ends. For each cut, one fragment is the more polar (the *head*)
//! and the other the less (the *tail*), and the cut scores
//!
//! ```text
//!     (mean polarity of head − mean polarity of tail) × atoms in tail
//! ```
//!
//! subject to the tail being at least two atoms, because a single atom is not a
//! phase. The molecule's amphiphilicity is the best score over all cuts.
//! Polarity is computed from the *parent's* bonds rather than the severed
//! fragments', so the bond joining head to tail still counts — which is right,
//! since in a real surfactant the linkage is part of what makes the head a head.
//!
//! # What is validated, and what is assumed
//!
//! The linear factor in tail size encodes **Traube's rule** (1891): each further
//! carbon strengthens a surfactant by roughly a constant increment, because the
//! hydrophobic driving force scales with the apolar surface that gets buried. As
//! with Graham's law in Phase 5b, be exact about the direction of the claim —
//! *the rule is input, not output*. It is put in deliberately, in one visible
//! factor, where it can be argued with.
//!
//! What the tests actually check is everything around it: that polarity is read
//! correctly off arbitrary discovered graphs, that molecules without polar
//! contrast score zero rather than noise, that rings are handled, that the
//! answer depends on the graph and not on the order its atoms happen to be
//! numbered in, and that the ranking of a network the engine *invented* matches
//! which molecules are actually surfactants.

use crate::element::PeriodicTable;
use crate::molecule::Molecule;

/// The polarity of every atom, in electronegativity units: for each atom, the
/// sum over its bonds of `order × |Δ electronegativity|`.
///
/// A nonpolar bond (carbon to carbon) contributes nothing; a strongly polar one
/// (oxygen to hydrogen) contributes to *both* of its atoms, which is the point —
/// polarity is a property of a region, not of an element.
pub fn atom_polarity(m: &Molecule, table: &PeriodicTable) -> Vec<f64> {
    let mut p = vec![0.0; m.atom_count()];
    for b in m.bonds() {
        let (ia, ib) = (b.a as usize, b.b as usize);
        let (za, zb) = (m.atoms()[ia].z, m.atoms()[ib].z);
        let (ea, eb) = match (table.get(za), table.get(zb)) {
            (Some(x), Some(y)) => (x, y),
            _ => continue,
        };
        let d = (ea.electronegativity() - eb.electronegativity()).abs() * b.order.slots() as f64;
        p[ia] += d;
        p[ib] += d;
    }
    p
}

/// The two atom sets a bond separates the molecule into, or `None` if removing
/// it leaves the molecule connected (a ring bond).
///
/// Deliberately parallel to [`crate::network::split`] — same flood fill, same
/// notion of a bridge — but returning *indices into the parent*, so the caller
/// can keep using the parent's bonds.
fn partition(m: &Molecule, bond_idx: usize) -> Option<(Vec<usize>, Vec<usize>)> {
    let n = m.atom_count();
    if bond_idx >= m.bonds().len() || n < 2 {
        return None;
    }
    let mut seen = vec![false; n];
    let mut stack = vec![0usize];
    seen[0] = true;
    while let Some(i) = stack.pop() {
        for (k, b) in m.bonds().iter().enumerate() {
            if k == bond_idx {
                continue;
            }
            let (a, c) = (b.a as usize, b.b as usize);
            let other = if a == i {
                c
            } else if c == i {
                a
            } else {
                continue;
            };
            if !seen[other] {
                seen[other] = true;
                stack.push(other);
            }
        }
    }
    let side_a: Vec<usize> = (0..n).filter(|&i| seen[i]).collect();
    let side_b: Vec<usize> = (0..n).filter(|&i| !seen[i]).collect();
    if side_b.is_empty() {
        None // still connected: the bond was part of a ring
    } else {
        Some((side_a, side_b))
    }
}

/// One way of cutting a molecule into a polar head and an apolar tail.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cut {
    /// Index of the bridge bond separating head from tail.
    pub bond: usize,
    /// Atoms in the more polar fragment.
    pub head_atoms: usize,
    /// Atoms in the less polar fragment — the part water would rather not see.
    pub tail_atoms: usize,
    /// Mean atom polarity of the head.
    pub head_density: f64,
    /// Mean atom polarity of the tail.
    pub tail_density: f64,
    /// `(head_density − tail_density) × tail_atoms`.
    pub score: f64,
}

/// The best head/tail cut of a molecule, or `None` if it has no usable bridge
/// bond at all — a lone atom, a pure ring, or nothing but single-atom tails.
pub fn best_cut(m: &Molecule, table: &PeriodicTable) -> Option<Cut> {
    let pol = atom_polarity(m, table);
    let mean = |ids: &[usize]| -> f64 {
        if ids.is_empty() {
            0.0
        } else {
            ids.iter().map(|&i| pol[i]).sum::<f64>() / ids.len() as f64
        }
    };

    let mut best: Option<Cut> = None;
    for k in 0..m.bonds().len() {
        let Some((sa, sb)) = partition(m, k) else { continue };
        let (da, db) = (mean(&sa), mean(&sb));
        // The head is whichever side is more polar. Ties break toward the
        // smaller fragment, so the answer never depends on where the flood fill
        // happened to start.
        let head_is_a = if da != db { da > db } else { sa.len() <= sb.len() };
        let (head, tail, dh, dt) =
            if head_is_a { (&sa, &sb, da, db) } else { (&sb, &sa, db, da) };

        // A single atom is not a phase. Snapping one terminal hydrogen off any
        // molecule leaves a slightly-less-polar remainder behind, so without
        // this every alkane in the world registers as a weak surfactant —
        // ethane scored 0.2 before a test caught it. A tail is the part that
        // gets *buried*, and one atom buries nothing. (Note a two-atom fragment
        // necessarily contains a heavy atom: hydrogen has one valence slot, so
        // a pair of hydrogens cannot both bond to each other and hang off the
        // rest of the molecule.)
        if tail.len() < 2 {
            continue;
        }
        let score = (dh - dt) * tail.len() as f64;
        let cut = Cut {
            bond: k,
            head_atoms: head.len(),
            tail_atoms: tail.len(),
            head_density: dh,
            tail_density: dt,
            score,
        };
        if best.map(|b| score > b.score).unwrap_or(true) {
            best = Some(cut);
        }
    }
    best
}

/// How strongly a molecule wants to sit at an interface, in
/// electronegativity-units × atoms. Zero for anything without polar contrast,
/// and zero — correctly — for a ring, which has no ends to be different.
///
/// This is a *shape* measure, not an energy. Turning it into a concentration at
/// which molecules actually assemble takes a declared scale, exactly as the
/// reaction barrier does; see [`crate::assembly`].
pub fn amphiphilicity(m: &Molecule, table: &PeriodicTable) -> f64 {
    best_cut(m, table).map(|c| c.score.max(0.0)).unwrap_or(0.0)
}
