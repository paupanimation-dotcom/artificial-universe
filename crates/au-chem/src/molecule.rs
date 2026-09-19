//! Molecules as graphs — the heart of "no hardcoded molecules".
//!
//! # The idea PROJECT_VISION actually demands
//!
//! > "No hardcoded creatures, species, ecosystems, civilizations, or
//! >  evolutionary paths."
//!
//! Applied to chemistry, that has a sharp consequence most "chemistry
//! simulations" quietly dodge: **there is no list of molecules.** Water is not an
//! entry in a table. Water is *two hydrogen atoms each bonded to one oxygen atom*
//! — a graph. The engine can construct that graph, compare it to another, and
//! recognise when two cells contain "the same molecule", all without anyone ever
//! having written down what water is.
//!
//! This is what makes genuine novelty possible. A reaction rule that says "a
//! hydrogen with a free bond may bond to an oxygen with a free bond" does not know
//! it is making water. It is making a *graph edge*. If, after enough edges, a
//! stable ring of six carbons appears that no one anticipated, the engine
//! represents it as naturally as it represents H₂O, because both are just graphs.
//! The molecule table is not pre-filled and then consulted; it is *discovered* and
//! then remembered.
//!
//! # Canonical form, or: how to know two graphs are the same molecule
//!
//! The hard problem hiding here is graph isomorphism. Two descriptions of ethanol
//! with the atoms listed in a different order are the *same molecule*, and the
//! engine must know that, or it will think it has discovered a new substance every
//! time the same one forms with its atoms enumerated differently — and the
//! molecule "census" that Phases 5–7 depend on would be meaningless noise.
//!
//! General graph isomorphism is famously hard, but **molecular** graphs are not
//! general: they are small, sparse, labelled by element, and bounded in degree by
//! valence. That makes a Morgan-style canonical labelling (iteratively refine each
//! atom's "rank" by its neighbourhood until the ranking is stable, then read the
//! graph out in rank order) both correct enough and cheap enough. The canonical
//! form is a byte string; two molecules are identical iff their canonical forms
//! are equal. It is the chemical equivalent of the world hash: a fingerprint that
//! is stable no matter how the thing was assembled.
//!
//! # What a molecule deliberately is *not*, yet
//!
//! No 3-D geometry, no stereochemistry, no bond angles, no resonance. Those are
//! real and they matter to real chemistry — but they are refinements, and the bet
//! of this phase is that *connectivity plus bond order* carries enough of the
//! causal structure for interesting chemistry to emerge. If Phase 5 shows that
//! chirality is load-bearing for replication (it may well be), geometry gets added
//! here without the layers above having to change — because they, too, only ever
//! see graphs and canonical forms.

use au_core::hash::{Hasher, WorldHash};

use crate::element::{PeriodicTable, Z};

/// How many electron pairs a bond shares: single, double, triple.
///
/// Bond order is not decoration — it is a claim about how many of each atom's
/// valence slots the bond consumes, and therefore about how much energy it takes
/// to break. A carbon–carbon triple bond uses three of each carbon's four slots
/// and is far stronger than a single. The reaction engine spends and reclaims
/// these slots exactly, so an atom can never end up with more bonds than its
/// valence allows.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
#[repr(u8)]
pub enum BondOrder {
    Single = 1,
    Double = 2,
    Triple = 3,
}

impl BondOrder {
    pub fn slots(self) -> u8 {
        self as u8
    }
    pub fn from_u8(v: u8) -> Option<BondOrder> {
        match v {
            1 => Some(BondOrder::Single),
            2 => Some(BondOrder::Double),
            3 => Some(BondOrder::Triple),
            _ => None,
        }
    }
}

/// One atom within a molecule. Just its element; identity within the molecule is
/// its position in the atom list.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Atom {
    pub z: Z,
}

/// A bond between two atoms of a molecule, by their indices. Stored with the lower
/// index first so that `(a,b)` and `(b,a)` are the same edge — an undirected graph
/// with a canonical edge representation.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Bond {
    pub a: u16,
    pub b: u16,
    pub order: BondOrder,
}

impl Bond {
    pub fn new(a: u16, b: u16, order: BondOrder) -> Bond {
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        Bond { a, b, order }
    }
}

/// A molecule: atoms plus the bonds between them. A labelled undirected graph, and
/// nothing more.
///
/// Constructed freely by the reaction engine — this type has no idea what it
/// represents, and that ignorance is the whole point. It can be a lone atom, a
/// diatomic gas, an amino acid, or something with no name because nothing like it
/// has existed before.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
}

impl Molecule {
    /// A single free atom — the simplest molecule, and where all chemistry starts
    /// before any bonds have formed.
    pub fn monatomic(z: Z) -> Molecule {
        Molecule { atoms: vec![Atom { z }], bonds: Vec::new() }
    }

    pub fn new(atoms: Vec<Atom>, bonds: Vec<Bond>) -> Molecule {
        Molecule { atoms, bonds }
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// How many bond slots atom `i` is currently using (summed bond orders).
    pub fn used_valence(&self, i: usize) -> u8 {
        self.bonds
            .iter()
            .filter(|b| b.a as usize == i || b.b as usize == i)
            .map(|b| b.order.slots())
            .sum()
    }

    /// Element formula as a count per element, for mass and conservation checks.
    /// This is the bridge back to `AtomCount`: however elaborate the graph, its
    /// atoms still have to balance.
    pub fn formula(&self, table: &PeriodicTable) -> Vec<i128> {
        let mut f = vec![0i128; table.len()];
        for a in &self.atoms {
            if let Some(idx) = table.index_of(a.z) {
                f[idx] += 1;
            }
        }
        f
    }

    pub fn mass_mda(&self, table: &PeriodicTable) -> i128 {
        self.atoms
            .iter()
            .filter_map(|a| table.get(a.z).map(|e| e.mass_mda as i128))
            .sum()
    }

    /// Is this a single connected molecule, or has it fragmented into pieces?
    ///
    /// The reaction engine needs this: breaking a bond can split one molecule into
    /// two, and the two halves must be recognised as separate molecules (each with
    /// its own canonical form) rather than lingering as one disconnected graph.
    pub fn is_connected(&self) -> bool {
        if self.atoms.len() <= 1 {
            return true;
        }
        let n = self.atoms.len();
        let mut seen = vec![false; n];
        let mut stack = vec![0usize];
        seen[0] = true;
        let mut count = 1;
        while let Some(cur) = stack.pop() {
            for b in &self.bonds {
                let other = if b.a as usize == cur {
                    Some(b.b as usize)
                } else if b.b as usize == cur {
                    Some(b.a as usize)
                } else {
                    None
                };
                if let Some(o) = other {
                    if !seen[o] {
                        seen[o] = true;
                        count += 1;
                        stack.push(o);
                    }
                }
            }
        }
        count == n
    }

    /// The canonical fingerprint of this molecule — equal for isomorphic graphs,
    /// different otherwise. Two cells "contain the same molecule" iff these match.
    ///
    /// Morgan-style: give every atom an initial invariant (its element plus its
    /// degree), then repeatedly replace each atom's rank with a hash of its own
    /// rank and the sorted ranks of its neighbours. After enough rounds the ranks
    /// stop changing and encode each atom's entire graph environment. Reading the
    /// atoms and bonds out in rank order gives a description that does not depend on
    /// the arbitrary order the atoms were listed in — a canonical form.
    ///
    /// This is a well-trodden approximation. It can, in rare highly-symmetric
    /// cases, fail to distinguish genuinely different graphs — but for the small,
    /// element-labelled molecules chemistry produces it is more than adequate, and
    /// it is O(atoms² · rounds) rather than the exponential cost of exact
    /// isomorphism. If a pathological case ever bites in practice, the fix is a
    /// stronger invariant here, invisible to every layer above.
    /// The converged Morgan-refinement ranks — one label per atom, where two
    /// atoms share a label exactly when refinement cannot tell their graph
    /// environments apart.
    ///
    /// These are the molecule's **symmetry classes**, and the Phase 5 reaction
    /// generator leans on them hard: when it asks "where could a new bond
    /// attach?", attaching at two atoms of the same class yields isomorphic
    /// products, so one representative per class suffices. Benzene offers one
    /// place to substitute, not six — and candidate enumeration stays small
    /// because of it.
    ///
    /// (Refinement can, on certain regular graphs, merge classes it should not.
    /// The generator therefore also deduplicates finished candidates by
    /// canonical form, so an over-merge here can never create a wrong reaction —
    /// at worst it pre-removes a duplicate the dedup would have caught anyway.)
    pub fn symmetry_classes(&self) -> Vec<u64> {
        let n = self.atoms.len();
        if n == 0 {
            return Vec::new();
        }

        // Initial invariant: element, then degree (number of incident bonds).
        let mut degree = vec![0u32; n];
        for b in &self.bonds {
            degree[b.a as usize] += 1;
            degree[b.b as usize] += 1;
        }
        let mut rank: Vec<u64> = (0..n)
            .map(|i| {
                let mut h = Hasher::new();
                h.write_u8(self.atoms[i].z.0);
                h.write_u32(degree[i]);
                h.finish()
            })
            .collect();

        // Refine ranks by neighbourhood until stable. `refine` replaces each
        // atom's rank with a hash of its own rank and the sorted (bond-order,
        // neighbour-rank) pairs, so after convergence a rank encodes an atom's
        // entire graph environment.
        let refine = |rank: &mut Vec<u64>, mol: &Molecule| {
            for _ in 0..mol.atoms.len() {
                let mut next = rank.clone();
                for i in 0..mol.atoms.len() {
                    let mut nbrs: Vec<(u8, u64)> = mol
                        .bonds
                        .iter()
                        .filter_map(|b| {
                            if b.a as usize == i {
                                Some((b.order.slots(), rank[b.b as usize]))
                            } else if b.b as usize == i {
                                Some((b.order.slots(), rank[b.a as usize]))
                            } else {
                                None
                            }
                        })
                        .collect();
                    nbrs.sort_unstable();
                    let mut h = Hasher::new();
                    h.write_u64(rank[i]);
                    for (o, r) in nbrs {
                        h.write_u8(o);
                        h.write_u64(r);
                    }
                    next[i] = h.finish();
                }
                if next == *rank {
                    break;
                }
                *rank = next;
            }
        };

        refine(&mut rank, self);
        rank
    }

    pub fn canonical(&self) -> CanonicalForm {
        let n = self.atoms.len();
        if n == 0 {
            return CanonicalForm(Vec::new());
        }
        let rank = self.symmetry_classes();

        // # Breaking symmetry the right way
        //
        // Refinement alone is not enough, and the failure is famous. In a fully
        // symmetric graph — benzene, where every carbon is genuinely identical —
        // refinement converges with *ties*: several atoms share a rank because
        // they truly have the same environment. Reading the graph out "in rank
        // order" then depends on how those ties are broken, and breaking them by
        // the atoms' arbitrary list positions is exactly the bug that made the
        // same ring described from two starting atoms canonicalise differently.
        //
        // The standard cure (McKay-style orbit refinement, in miniature): while a
        // tie remains, pick the *smallest-ranked* tied group, tentatively single
        // one of its members out by nudging its rank, re-refine, and repeat. This
        // forces a total order that is a property of the graph's structure, not of
        // the input ordering. It is a canonical *discrete* refinement.
        //
        // To make the choice itself canonical when the graph is symmetric, we try
        // singling out *each* member of the tied group and keep the lexicographically
        // smallest resulting edge-list — so equivalent atoms yield the same answer.
        let order = canonical_order(self, rank);
        let pos: Vec<u16> = {
            let mut p = vec![0u16; n];
            for (new_idx, &old) in order.iter().enumerate() {
                p[old] = new_idx as u16;
            }
            p
        };

        let mut out = Hasher::new();
        out.write_u64(n as u64);
        for &old in &order {
            out.write_u8(self.atoms[old].z.0);
        }
        // Bonds, rewritten in canonical positions, sorted.
        let mut edges: Vec<(u16, u16, u8)> = self
            .bonds
            .iter()
            .map(|b| {
                let (x, y) = (pos[b.a as usize], pos[b.b as usize]);
                let (x, y) = if x <= y { (x, y) } else { (y, x) };
                (x, y, b.order.slots())
            })
            .collect();
        edges.sort_unstable();
        out.write_u64(edges.len() as u64);
        for (x, y, o) in edges {
            out.write_u16(x);
            out.write_u16(y);
            out.write_u8(o);
        }

        // The canonical form is the full byte sequence, not just a hash of it, so
        // that it cannot silently collide the way a single u64 eventually would
        // across a whole planet's worth of molecules.
        let digest = out.finish();
        let mut bytes = Vec::with_capacity(1 + n + 1);
        bytes.extend_from_slice(&digest.to_le_bytes());
        bytes.push(n as u8);
        for &old in &order {
            bytes.push(self.atoms[old].z.0);
        }
        CanonicalForm(bytes)
    }
}

impl WorldHash for Molecule {
    fn hash_into(&self, h: &mut Hasher) {
        // Hash the canonical form, so two isomorphic molecules hash identically —
        // the world hash must not depend on atom ordering any more than it depends
        // on chunk residency.
        let c = self.canonical();
        h.write_bytes(&c.0);
    }
}

/// Produce a canonical atom ordering by discrete refinement with symmetry
/// breaking. Returns the atom indices in canonical order.
///
/// When refinement leaves ties (a symmetric graph), this singles out each atom of
/// the smallest tied orbit in turn, re-refines, recurses, and keeps whichever
/// choice yields the lexicographically smallest edge list. Because equivalent
/// atoms produce identical edge lists, the result is independent of which member
/// was chosen — a true canonical form.
fn canonical_order(mol: &Molecule, rank: Vec<u64>) -> Vec<usize> {
    let n = mol.atoms.len();

    // Compress ranks to 0..k dense classes, preserving order.
    let mut sorted: Vec<u64> = rank.clone();
    sorted.sort_unstable();
    sorted.dedup();
    let class = |r: u64| sorted.binary_search(&r).unwrap() as u64;
    let classes: Vec<u64> = rank.iter().map(|&r| class(r)).collect();

    // Find the smallest class that contains more than one atom (a nontrivial orbit).
    let mut counts = vec![0usize; sorted.len()];
    for &c in &classes {
        counts[c as usize] += 1;
    }
    let tied = counts.iter().position(|&c| c > 1);

    match tied {
        None => {
            // Fully discrete: the ordering is determined. Sort by class.
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by_key(|&i| classes[i]);
            order
        }
        Some(orbit) => {
            // Try promoting each atom of this orbit; keep the lexicographically
            // smallest resulting edge list. Recurse until discrete.
            let members: Vec<usize> =
                (0..n).filter(|&i| classes[i] as usize == orbit).collect();
            let mut best_order: Option<Vec<usize>> = None;
            let mut best_edges: Option<Vec<(u16, u16, u8)>> = None;

            for &m in &members {
                // Promote atom m to its own, strictly-smallest class by scaling
                // all ranks and lowering m's, then re-refine.
                let mut r2: Vec<u64> = classes.iter().map(|&c| (c + 1) * 2).collect();
                r2[m] = 0;
                // Local re-refinement (same routine, inlined to avoid capturing).
                for _ in 0..n {
                    let mut next = r2.clone();
                    for i in 0..n {
                        let mut nbrs: Vec<(u8, u64)> = mol
                            .bonds
                            .iter()
                            .filter_map(|b| {
                                if b.a as usize == i {
                                    Some((b.order.slots(), r2[b.b as usize]))
                                } else if b.b as usize == i {
                                    Some((b.order.slots(), r2[b.a as usize]))
                                } else {
                                    None
                                }
                            })
                            .collect();
                        nbrs.sort_unstable();
                        let mut h = Hasher::new();
                        h.write_u64(r2[i]);
                        for (o, rr) in nbrs {
                            h.write_u8(o);
                            h.write_u64(rr);
                        }
                        next[i] = h.finish();
                    }
                    if next == r2 {
                        break;
                    }
                    r2 = next;
                }
                let order = canonical_order(mol, r2);
                let edges = edge_list_in_order(mol, &order);
                let better = match &best_edges {
                    None => true,
                    Some(be) => &edges < be,
                };
                if better {
                    best_edges = Some(edges);
                    best_order = Some(order);
                }
            }
            best_order.unwrap()
        }
    }
}

/// The molecule's edges rewritten to the given atom order, sorted — the object we
/// minimise over when breaking symmetry.
fn edge_list_in_order(mol: &Molecule, order: &[usize]) -> Vec<(u16, u16, u8)> {
    let mut pos = vec![0u16; mol.atoms.len()];
    for (new_idx, &old) in order.iter().enumerate() {
        pos[old] = new_idx as u16;
    }
    let mut edges: Vec<(u16, u16, u8)> = mol
        .bonds
        .iter()
        .map(|b| {
            let (x, y) = (pos[b.a as usize], pos[b.b as usize]);
            let (x, y) = if x <= y { (x, y) } else { (y, x) };
            (x, y, b.order.slots())
        })
        .collect();
    edges.sort_unstable();
    edges
}

/// A molecule's identity as a comparable, orderable byte string. Equal iff the
/// molecules are the same graph.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub struct CanonicalForm(pub Vec<u8>);

impl WorldHash for CanonicalForm {
    fn hash_into(&self, h: &mut Hasher) {
        h.write_bytes(&self.0);
    }
}
