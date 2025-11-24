pub(crate) type Float = f64;

pub const DELTA_T: Float = 1.0;

pub const DELTA_X: Float = 1.0;

pub const LATTICE_DENSITY: Float = 1.0;

pub const CS_2: Float = 1.0 / 3.0 * DELTA_X * DELTA_X / DELTA_T / DELTA_T;

pub const CS_2_INV: Float = 3.0;

pub const CS_4_INV: Float = 9.0;

pub(crate) const MIN_ITER: usize = 10;

pub(crate) const TOLERANCE_DENSITY: Float = 1e-7;

pub(crate) const TOLERANCE_VELOCITY: Float = 1e-7;

pub(crate) const TOLERANCE_SCALAR_VALUE: Float = 1e-7;
