//! The layer stack, encoded as a type the engine can actually check.
//!
//! ARCHITECTURE states the dependency chain and the Golden Rule:
//!
//! > "No higher-level system should directly create results that should come
//! >  from lower-level systems."
//!
//! In most codebases that sentence lives in a document, and the document
//! quietly becomes a lie over about eighteen months. Here it is a boot-time
//! assertion.
//!
//! Every system declares which layer it *writes* and which layers it *reads*.
//! The scheduler rejects, loudly and at startup, any system that reads a layer
//! above its own. You cannot write an ecology system that reaches up into
//! civilization for an answer, because the engine will not start.
//!
//! The ordering of this enum *is* the architecture. Adding a variant in the
//! wrong place changes what the universe is allowed to do.

/// The layers of the simulation, in strict dependency order.
///
/// Lower discriminant = more fundamental. Derived `Ord` gives us the check.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
#[repr(u8)]
pub enum Layer {
    /// Constants, seeds, celestial structure. Nothing is below this.
    Universe = 0,
    Physics = 1,
    Matter = 2,
    Chemistry = 3,
    Planet = 4,
    Environment = 5,
    /// The chemistry→life transition. Deliberately its own layer: it is the
    /// hardest boundary in the project and it must not be allowed to hide
    /// inside either neighbour.
    LifeChemistry = 6,
    Genetics = 7,
    Evolution = 8,
    Organisms = 9,
    Behavior = 10,
    Ecology = 11,
    Intelligence = 12,
    Civilization = 13,
}

impl Layer {
    pub const ALL: [Layer; 14] = [
        Layer::Universe,
        Layer::Physics,
        Layer::Matter,
        Layer::Chemistry,
        Layer::Planet,
        Layer::Environment,
        Layer::LifeChemistry,
        Layer::Genetics,
        Layer::Evolution,
        Layer::Organisms,
        Layer::Behavior,
        Layer::Ecology,
        Layer::Intelligence,
        Layer::Civilization,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Layer::Universe => "Universe",
            Layer::Physics => "Physics",
            Layer::Matter => "Matter",
            Layer::Chemistry => "Chemistry",
            Layer::Planet => "Planet",
            Layer::Environment => "Environment",
            Layer::LifeChemistry => "LifeChemistry",
            Layer::Genetics => "Genetics",
            Layer::Evolution => "Evolution",
            Layer::Organisms => "Organisms",
            Layer::Behavior => "Behavior",
            Layer::Ecology => "Ecology",
            Layer::Intelligence => "Intelligence",
            Layer::Civilization => "Civilization",
        }
    }

    /// May a system that writes `self` legally read `other`?
    ///
    /// Yes if and only if `other` is at or below `self`. A system may read its
    /// own layer (organisms see other organisms) and everything beneath it
    /// (organisms feel gravity). It may never read upward: that would be the
    /// engine asking "what content do I need?" instead of "what happens?".
    #[inline]
    pub fn may_read(self, other: Layer) -> bool {
        other <= self
    }
}
