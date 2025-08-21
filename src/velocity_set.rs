// ------------------------------------------------------------------------------- MODULES

pub mod d2q9;
pub mod d3q15;
pub mod d3q19;
pub mod d3q27;

// ------------------------------------------------------------------------------- IMPORTS

use crate::prelude::*;
use crate::{FACES_2D, FACES_3D};
use std::collections::HashMap;

pub type VelocityComputation = fn(Float, Vec<Float>) -> Vec<Float>;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VelocitySet {
    D2Q9,
    D3Q15,
    D3Q19,
    D3Q27,
}

impl VelocitySet {
    pub fn get_velocity_set_parameters(&self) -> VelocitySetParameters {
        match self {
            D2Q9 => VelocitySetParameters {
                velocity_set: D2Q9,
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
                    (West, 1),
                    (East, 3),
                    (South, 2),
                    (North, 4),
                ]),
                velocity_computation: Some(d2q9::velocity_computation),
            },
            D3Q15 => VelocitySetParameters {
                velocity_set: D3Q15,
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
                    (West, 1),
                    (East, 2),
                    (South, 3),
                    (North, 4),
                    (Bottom, 5),
                    (Top, 6),
                ]),
                velocity_computation: None,
            },
            D3Q19 => VelocitySetParameters {
                velocity_set: D3Q19,
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
                    (West, 1),
                    (East, 2),
                    (South, 3),
                    (North, 4),
                    (Bottom, 5),
                    (Top, 6),
                ]),
                velocity_computation: None,
            },
            D3Q27 => VelocitySetParameters {
                velocity_set: D3Q27,
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
                    (West, 1),
                    (East, 2),
                    (South, 3),
                    (North, 4),
                    (Bottom, 5),
                    (Top, 6),
                ]),
                velocity_computation: Some(d3q27::velocity_computation),
            },
        }
    }
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
    pub velocity_computation: Option<VelocityComputation>,
}

impl Default for VelocitySetParameters {
    fn default() -> Self {
        D2Q9.get_velocity_set_parameters()
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
            3 => D3Q27.get_velocity_set_parameters(),
            _ => panic!("Unsupported dimension: {dim}"),
        }
    }
}
