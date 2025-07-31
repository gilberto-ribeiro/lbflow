pub mod constants;
pub mod io;
pub mod momentum;
mod parameters;
mod simulation;
mod thermal;
pub mod velocity_set;

use velocity_set::*;

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum BoundaryFace {
    West = 0,
    East = 1,
    South = 2,
    North = 3,
    Bottom = 4,
    Top = 5,
}

const FACES_2D: [BoundaryFace; 4] = [
    BoundaryFace::West,
    BoundaryFace::East,
    BoundaryFace::South,
    BoundaryFace::North,
];

const FACES_3D: [BoundaryFace; 6] = [
    BoundaryFace::West,
    BoundaryFace::East,
    BoundaryFace::South,
    BoundaryFace::North,
    BoundaryFace::Bottom,
    BoundaryFace::Top,
];

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum NodeType {
    Fluid = 0,
    Solid = 1,
}
