//! Where physics meets storage.
//!
//! `au-physics` knows nothing about chunks, and `au-data` knows nothing about
//! energy. This file is the only place that knows both, and it stays small on
//! purpose — the moment it starts making physical decisions, the layering has
//! sprung a leak.
//!
//! # Column ids are permanent
//!
//! A save file names its columns by number. Renumbering one does not fail to
//! load — it loads *wrong*, silently reinterpreting a mass as an energy. So the
//! physics layer owns the `0x1xxx` range and these three numbers are now
//! immortal. Chemistry gets `0x2xxx`, genetics `0x3xxx`, and so on.

use au_data::column::{ColumnId, ColumnRegistry, ColumnSet, VecColumn};

/// Mass per cell, micrograms. Exact. Constant in Phase 2 — nothing moves matter
/// yet — and the day advection arrives it gets the same symmetric-flux treatment
/// energy already has, for the same reason.
pub const PHYS_MASS: ColumnId = ColumnId(0x1001);

/// Energy per cell, microjoules. Exact. The only thing Phase 2 changes.
pub const PHYS_ENERGY: ColumnId = ColumnId(0x1002);

/// Which substance. Index into the world's `MaterialRegistry`.
pub const PHYS_MATERIAL: ColumnId = ColumnId(0x1003);

/// Momentum on the −x / −y / −z face of each cell. Picogram-metres per second.
///
/// On *faces*, not at cell centres. Co-locating velocity and pressure decouples
/// the odd and even cells in the pressure equation and produces an invisible
/// checkerboard that drives fictitious flow — see `au_physics::grid::FluidView`.
pub const PHYS_MOM_X: ColumnId = ColumnId(0x1004);
pub const PHYS_MOM_Y: ColumnId = ColumnId(0x1005);
pub const PHYS_MOM_Z: ColumnId = ColumnId(0x1006);

/// Dynamic pressure. Solver state, not a fundamental field — but it *is* state,
/// because the fixed-iteration solve is warm-started from it, and a resumed world
/// must get the same warm start a continuous one would have had.
pub const PHYS_PRESSURE: ColumnId = ColumnId(0x1007);

pub fn register(reg: &mut ColumnRegistry) {
    reg.register::<i128>(PHYS_MASS);
    reg.register::<i128>(PHYS_ENERGY);
    reg.register::<u16>(PHYS_MATERIAL);
    reg.register::<i128>(PHYS_MOM_X);
    reg.register::<i128>(PHYS_MOM_Y);
    reg.register::<i128>(PHYS_MOM_Z);
    reg.register::<f64>(PHYS_PRESSURE);
}

/// Does this chunk carry a velocity field?
pub fn has_fluid(cols: &ColumnSet) -> bool {
    cols.get::<i128>(PHYS_MOM_X).is_some() && cols.get::<f64>(PHYS_PRESSURE).is_some()
}

/// Give a chunk a velocity field. At rest.
pub fn install_fluid(cols: &mut ColumnSet, cells: usize) {
    for id in [PHYS_MOM_X, PHYS_MOM_Y, PHYS_MOM_Z] {
        cols.insert(Box::new(VecColumn::<i128>::with_data(id, vec![0; cells])));
    }
    cols.insert(Box::new(VecColumn::<f64>::with_data(PHYS_PRESSURE, vec![0.0; cells])));
}

/// Does this chunk have physics in it?
///
/// Most chunks will not. A planet is mostly vacuum and the parts of it worth
/// integrating in detail are a thin skin.
pub fn has_physics(cols: &ColumnSet) -> bool {
    cols.get::<i128>(PHYS_MASS).is_some()
        && cols.get::<i128>(PHYS_ENERGY).is_some()
        && cols.get::<u16>(PHYS_MATERIAL).is_some()
}

/// Install empty physics columns of the right length.
///
/// Called by generators, never by systems: a system that could *create* matter
/// where there was none is a system that can conjure, and conjuring is the thing
/// the Golden Rule exists to prevent.
pub fn install(cols: &mut ColumnSet, cells: usize) {
    cols.insert(Box::new(VecColumn::<i128>::with_data(PHYS_MASS, vec![0; cells])));
    cols.insert(Box::new(VecColumn::<i128>::with_data(PHYS_ENERGY, vec![0; cells])));
    cols.insert(Box::new(VecColumn::<u16>::with_data(
        PHYS_MATERIAL,
        vec![au_physics::MaterialId::VACUUM.0; cells],
    )));
}
