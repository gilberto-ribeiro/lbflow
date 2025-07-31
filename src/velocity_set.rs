pub mod d2q9;
pub mod d3q15;
pub mod d3q19;
pub mod d3q27;

use crate::constants::Float;
use crate::{BoundaryFace, FACES_2D, FACES_3D};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VelocitySet {
    D2Q9,
    D3Q15,
    D3Q19,
    D3Q27,
}

#[derive(Debug, PartialEq)]
pub struct VelocitySetParameters {
    pub velocity_set: VelocitySet,
    pub d: usize,
    pub q: usize,
    pub c: Vec<Vec<i32>>,
    pub w: Vec<Float>,
    pub q_bar: Vec<usize>,
    pub q_faces: HashMap<BoundaryFace, Vec<usize>>,
    pub face_normal_directions: HashMap<BoundaryFace, usize>,
    pub velocity_computation: Option<fn(Float, Vec<Float>) -> Vec<Float>>,
}

impl Default for VelocitySetParameters {
    fn default() -> Self {
        get_velocity_set_parameters(VelocitySet::D2Q9)
    }
}

/// # Examples
/// ## 2D:
/// ```
/// # use lbflow::velocity_set::VelocitySet;
/// # use lbflow::velocity_set::VelocitySetParameters;
/// let velocity_set_parameters = VelocitySetParameters::test_default(2);
/// assert_eq!(velocity_set_parameters.d, 2);
/// assert_eq!(velocity_set_parameters.q, 9);
/// assert_eq!(velocity_set_parameters.c.len(), 9);
/// assert_eq!(velocity_set_parameters.w.len(), 9);
/// ```
/// ## 3D:
/// ```
/// # use lbflow::velocity_set::VelocitySet;
/// # use lbflow::velocity_set::VelocitySetParameters;
/// let velocity_set_parameters = VelocitySetParameters::test_default(3);
/// assert_eq!(velocity_set_parameters.d, 3);
/// assert_eq!(velocity_set_parameters.q, 27);
/// assert_eq!(velocity_set_parameters.c.len(), 27);
/// assert_eq!(velocity_set_parameters.w.len(), 27);
/// ```
impl VelocitySetParameters {
    pub fn test_default(dim: usize) -> Self {
        match dim {
            2 => Default::default(),
            3 => get_velocity_set_parameters(VelocitySet::D3Q27),
            _ => panic!("Unsupported dimension: {dim}"),
        }
    }
}

pub fn get_velocity_set_parameters(velocity_set: VelocitySet) -> VelocitySetParameters {
    match velocity_set {
        VelocitySet::D2Q9 => VelocitySetParameters {
            velocity_set: VelocitySet::D2Q9,
            d: d2q9::D,
            q: d2q9::Q,
            c: d2q9::C.iter().map(|&arr| arr.to_vec()).collect(),
            w: d2q9::W.to_vec(),
            q_bar: d2q9::Q_BAR.to_vec(),
            q_faces: FACES_2D
                .iter()
                .zip(d2q9::Q_FACES.iter())
                .map(|(face, &arr)| (*face, arr.to_vec()))
                .collect(),
            face_normal_directions: HashMap::from([
                (BoundaryFace::West, 1),
                (BoundaryFace::East, 3),
                (BoundaryFace::South, 2),
                (BoundaryFace::North, 4),
            ]),
            velocity_computation: Some(d2q9::velocity_computation),
        },
        VelocitySet::D3Q15 => VelocitySetParameters {
            velocity_set: VelocitySet::D3Q15,
            d: d3q15::D,
            q: d3q15::Q,
            c: d3q15::C.iter().map(|&arr| arr.to_vec()).collect(),
            w: d3q15::W.to_vec(),
            q_bar: d3q15::Q_BAR.to_vec(),
            q_faces: FACES_3D
                .iter()
                .zip(d3q15::Q_FACES.iter())
                .map(|(face, &arr)| (*face, arr.to_vec()))
                .collect(),
            face_normal_directions: HashMap::from([
                (BoundaryFace::West, 1),
                (BoundaryFace::East, 2),
                (BoundaryFace::South, 3),
                (BoundaryFace::North, 4),
                (BoundaryFace::Bottom, 5),
                (BoundaryFace::Top, 6),
            ]),
            velocity_computation: None,
        },
        VelocitySet::D3Q19 => VelocitySetParameters {
            velocity_set: VelocitySet::D3Q19,
            d: d3q19::D,
            q: d3q19::Q,
            c: d3q19::C.iter().map(|&arr| arr.to_vec()).collect(),
            w: d3q19::W.to_vec(),
            q_bar: d3q19::Q_BAR.to_vec(),
            q_faces: FACES_3D
                .iter()
                .zip(d3q19::Q_FACES.iter())
                .map(|(face, &arr)| (*face, arr.to_vec()))
                .collect(),
            face_normal_directions: HashMap::from([
                (BoundaryFace::West, 1),
                (BoundaryFace::East, 2),
                (BoundaryFace::South, 3),
                (BoundaryFace::North, 4),
                (BoundaryFace::Bottom, 5),
                (BoundaryFace::Top, 6),
            ]),
            velocity_computation: None,
        },
        VelocitySet::D3Q27 => VelocitySetParameters {
            velocity_set: VelocitySet::D3Q27,
            d: d3q27::D,
            q: d3q27::Q,
            c: d3q27::C.iter().map(|&arr| arr.to_vec()).collect(),
            w: d3q27::W.to_vec(),
            q_bar: d3q27::Q_BAR.to_vec(),
            q_faces: FACES_3D
                .iter()
                .zip(d3q27::Q_FACES.iter())
                .map(|(face, &arr)| (*face, arr.to_vec()))
                .collect(),
            face_normal_directions: HashMap::from([
                (BoundaryFace::West, 1),
                (BoundaryFace::East, 2),
                (BoundaryFace::South, 3),
                (BoundaryFace::North, 4),
                (BoundaryFace::Bottom, 5),
                (BoundaryFace::Top, 6),
            ]),
            velocity_computation: Some(d3q27::velocity_computation),
        },
    }
}
