// ------------------------------------------------------------------------------- MODULES

pub mod d2q9;
pub mod d3q15;
pub mod d3q19;
pub mod d3q27;

// ------------------------------------------------------------------------------- IMPORTS

use crate::prelude::*;
use crate::{FACES_2D, FACES_3D};
use std::collections::HashMap;

pub type VectorComputation = fn(Float, Vec<Float>) -> Vec<Float>;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VelocitySet {
    D2Q9,
    D3Q15,
    D3Q19,
    D3Q27,
}

impl VelocitySet {
    pub fn get_velocity_set_parameters(&self) -> Parameters {
        match self {
            D2Q9 => Parameters {
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
                mrt_matrix: d2q9::MRT_MATRIX.iter().map(|row| row.to_vec()).collect(),
                mrt_inverse_matrix: d2q9::MRT_INVERSE_MATRIX
                    .iter()
                    .map(|row| row.to_vec())
                    .collect(),
                mrt_equilibrium_moments_computation: Some(
                    d2q9::mrt_equilibrium_moments_computation,
                ),
            },
            D3Q15 => Parameters {
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
                mrt_matrix: d3q15::MRT_MATRIX.iter().map(|row| row.to_vec()).collect(),
                mrt_inverse_matrix: d3q15::MRT_INVERSE_MATRIX
                    .iter()
                    .map(|row| row.to_vec())
                    .collect(),
                mrt_equilibrium_moments_computation: Some(
                    d3q15::mrt_equilibrium_moments_computation,
                ),
            },
            D3Q19 => Parameters {
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
                mrt_matrix: d3q19::MRT_MATRIX.iter().map(|row| row.to_vec()).collect(),
                mrt_inverse_matrix: d3q19::MRT_INVERSE_MATRIX
                    .iter()
                    .map(|row| row.to_vec())
                    .collect(),
                mrt_equilibrium_moments_computation: Some(
                    d3q19::mrt_equilibrium_moments_computation,
                ),
            },
            D3Q27 => Parameters {
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
                mrt_matrix: d3q27::MRT_MATRIX.iter().map(|row| row.to_vec()).collect(),
                mrt_inverse_matrix: compute_inverse_matrix(
                    d3q27::MRT_MATRIX.iter().map(|row| row.to_vec()).collect(),
                ),
                mrt_equilibrium_moments_computation: None,
            },
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Parameters {
    pub velocity_set: VelocitySet,
    pub d: usize,
    pub q: usize,
    pub c: Vec<Vec<i32>>,
    pub w: Vec<Float>,
    pub q_bar: Vec<usize>,
    pub q_faces: HashMap<BoundaryFace, Vec<usize>>,
    pub face_normal_directions: HashMap<BoundaryFace, usize>,
    pub velocity_computation: Option<VectorComputation>,
    pub mrt_matrix: Vec<Vec<Float>>,
    pub mrt_inverse_matrix: Vec<Vec<Float>>,
    pub mrt_equilibrium_moments_computation: Option<VectorComputation>,
}

impl Default for Parameters {
    fn default() -> Self {
        D2Q9.get_velocity_set_parameters()
    }
}

/// # Examples
/// ## 2D:
/// ```
/// # use lbflow::velocity_set::VelocitySet;
/// # use lbflow::velocity_set::Parameters;
/// let velocity_set_parameters = Parameters::test_default(2);
/// assert_eq!(velocity_set_parameters.d, 2);
/// assert_eq!(velocity_set_parameters.q, 9);
/// assert_eq!(velocity_set_parameters.c.len(), 9);
/// assert_eq!(velocity_set_parameters.w.len(), 9);
/// ```
/// ## 3D:
/// ```
/// # use lbflow::velocity_set::VelocitySet;
/// # use lbflow::velocity_set::Parameters;
/// let velocity_set_parameters = Parameters::test_default(3);
/// assert_eq!(velocity_set_parameters.d, 3);
/// assert_eq!(velocity_set_parameters.q, 27);
/// assert_eq!(velocity_set_parameters.c.len(), 27);
/// assert_eq!(velocity_set_parameters.w.len(), 27);
/// ```
impl Parameters {
    pub fn test_default(dim: usize) -> Self {
        match dim {
            2 => Default::default(),
            3 => D3Q27.get_velocity_set_parameters(),
            _ => panic!("Unsupported dimension: {dim}"),
        }
    }
}

impl Parameters {
    pub fn get_d(&self) -> usize {
        self.d
    }

    pub fn get_q(&self) -> usize {
        self.q
    }

    pub fn get_c(&self) -> &Vec<Vec<i32>> {
        &self.c
    }

    pub fn get_w(&self) -> &Vec<Float> {
        &self.w
    }

    pub fn get_q_bar(&self) -> &Vec<usize> {
        &self.q_bar
    }

    pub fn get_q_faces(&self, boundary_face: &BoundaryFace) -> &Vec<usize> {
        self.q_faces
            .get(boundary_face)
            .expect("Boundary face not found")
    }

    pub fn get_face_normal_direction(&self, boundary_face: &BoundaryFace) -> &usize {
        self.face_normal_directions
            .get(boundary_face)
            .expect("Boundary face not found")
    }

    pub fn get_velocity_computation(&self) -> Option<VectorComputation> {
        self.velocity_computation
    }

    pub fn get_opposite_direction(&self, direction: usize) -> usize {
        self.get_q_bar()[direction]
    }

    pub fn get_mrt_matrix(&self) -> &Vec<Vec<Float>> {
        &self.mrt_matrix
    }

    pub fn get_mrt_inverse_matrix(&self) -> &Vec<Vec<Float>> {
        &self.mrt_inverse_matrix
    }

    pub fn get_mrt_equilibrium_moments_computation(&self) -> Option<&VectorComputation> {
        self.mrt_equilibrium_moments_computation.as_ref()
    }
}

fn compute_inverse_matrix(matrix: Vec<Vec<Float>>) -> Vec<Vec<Float>> {
    let rows = matrix.len();
    let cols = matrix[0].len();
    let flat_data: Vec<f64> = matrix.iter().flatten().cloned().collect();
    let dmatrix = nalgebra::DMatrix::from_row_slice(rows, cols, &flat_data);
    match dmatrix.try_inverse() {
        Some(inverse) => (0..inverse.nrows())
            .map(|i| inverse.row(i).iter().cloned().collect())
            .collect::<Vec<Vec<f64>>>(),
        None => {
            panic!("Singular matrix!")
        }
    }
}
