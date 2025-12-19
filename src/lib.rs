// ------------------------------------------------------------------------------- MODULES

mod cli;
mod constants;
mod io;
mod kernel;
pub mod momentum;
pub mod passive_scalar;
pub mod prelude;
mod prelude_crate;
mod velocity_set;

// ------------------------------------------------------------------------------- IMPORTS

use prelude_crate::*;

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum BoundaryFace {
    West = 0,
    East = 1,
    South = 2,
    North = 3,
    Bottom = 4,
    Top = 5,
}

const FACES_2D: [BoundaryFace; 4] = [West, East, South, North];

const FACES_3D: [BoundaryFace; 6] = [West, East, South, North, Bottom, Top];

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum NodeType {
    Fluid = 0,
    Solid = 1,
}

#[derive(Debug)]
pub enum CollisionOperator {
    BGK(Float),
    TRT(Float, Float),
    MRT(Vec<Float>),
}

impl Default for CollisionOperator {
    fn default() -> Self {
        BGK(0.5)
    }
}
