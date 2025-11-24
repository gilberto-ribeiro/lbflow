// ------------------------------------------------------------------------------- MODULES

mod d2q9;
mod d3q15;
mod d3q19;
mod d3q27;

// ------------------------------------------------------------------------------- IMPORTS

use crate::prelude_crate::*;
use crate::{FACES_2D, FACES_3D};
use std::collections::HashMap;

pub(crate) type VectorComputation = fn(Float, Vec<Float>) -> Vec<Float>;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VelocitySet {
    D2Q9 = 0,
    D3Q15 = 1,
    D3Q19 = 2,
    D3Q27 = 3,
}

impl VelocitySet {
    pub(crate) fn get_velocity_set_parameters(&self) -> Parameters {
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
pub(crate) struct Parameters {
    pub(crate) velocity_set: VelocitySet,
    pub(crate) d: usize,
    pub(crate) q: usize,
    pub(crate) c: Vec<Vec<i32>>,
    pub(crate) w: Vec<Float>,
    pub(crate) q_bar: Vec<usize>,
    pub(crate) q_faces: HashMap<BoundaryFace, Vec<usize>>,
    pub(crate) face_normal_directions: HashMap<BoundaryFace, usize>,
    pub(crate) velocity_computation: Option<VectorComputation>,
    pub(crate) mrt_matrix: Vec<Vec<Float>>,
    pub(crate) mrt_inverse_matrix: Vec<Vec<Float>>,
    pub(crate) mrt_equilibrium_moments_computation: Option<VectorComputation>,
}

impl Default for Parameters {
    fn default() -> Self {
        D2Q9.get_velocity_set_parameters()
    }
}

impl Parameters {
    pub(crate) fn test_default(dim: usize) -> Self {
        match dim {
            2 => Default::default(),
            3 => D3Q19.get_velocity_set_parameters(),
            _ => panic!("Unsupported dimension: {dim}"),
        }
    }
}

impl Parameters {
    pub(crate) fn get_d(&self) -> usize {
        self.d
    }

    pub(crate) fn get_q(&self) -> usize {
        self.q
    }

    pub(crate) fn get_c(&self) -> &Vec<Vec<i32>> {
        &self.c
    }

    pub(crate) fn get_w(&self) -> &Vec<Float> {
        &self.w
    }

    pub(crate) fn get_q_bar(&self) -> &Vec<usize> {
        &self.q_bar
    }

    pub(crate) fn get_q_faces(&self, boundary_face: &BoundaryFace) -> &Vec<usize> {
        self.q_faces
            .get(boundary_face)
            .expect("Boundary face not found")
    }

    pub(crate) fn get_face_normal_direction(&self, boundary_face: &BoundaryFace) -> usize {
        *self
            .face_normal_directions
            .get(boundary_face)
            .expect("Boundary face not found")
    }

    pub(crate) fn get_velocity_computation(&self) -> Option<VectorComputation> {
        self.velocity_computation
    }

    pub(crate) fn get_opposite_direction(&self, direction: usize) -> usize {
        self.get_q_bar()[direction]
    }

    pub(crate) fn get_mrt_matrix(&self) -> &Vec<Vec<Float>> {
        &self.mrt_matrix
    }

    pub(crate) fn get_mrt_inverse_matrix(&self) -> &Vec<Vec<Float>> {
        &self.mrt_inverse_matrix
    }

    pub(crate) fn get_mrt_equilibrium_moments_computation(&self) -> Option<&VectorComputation> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_d_d2q9() {
        let vel_set_params = Parameters::test_default(2);

        assert_eq!(vel_set_params.get_d(), 2);
    }

    #[test]
    fn test_get_d_d3q19() {
        let vel_set_params = Parameters::test_default(3);

        assert_eq!(vel_set_params.get_d(), 3);
    }

    #[test]
    fn test_get_q_d2q9() {
        let vel_set_params = Parameters::test_default(2);

        assert_eq!(vel_set_params.get_q(), 9);
    }

    #[test]
    fn test_get_q_d3q19() {
        let vel_set_params = Parameters::test_default(3);

        assert_eq!(vel_set_params.get_q(), 19);
    }

    #[test]
    fn test_get_c_d2q9() {
        let vel_set_params = Parameters::test_default(2);

        let c = vel_set_params.get_c();

        assert_eq!(c[0], vec![0, 0]);
        assert_eq!(c[1], vec![1, 0]);
        assert_eq!(c[5], vec![1, 1]);
    }

    #[test]
    fn test_get_c_d3q19() {
        let vel_set_params = Parameters::test_default(3);

        let c = vel_set_params.get_c();

        assert_eq!(c[0], vec![0, 0, 0]);
        assert_eq!(c[1], vec![1, 0, 0]);
        assert_eq!(c[7], vec![1, 1, 0]);
        assert_eq!(c[9], vec![1, 0, 1]);
    }

    #[test]
    fn test_get_w_d2q9() {
        let vel_set_params = Parameters::test_default(2);

        let w = vel_set_params.get_w();

        assert!(w[0] - 4.0 / 9.0 < 1e-12);
        assert!(w[1] - 1.0 / 9.0 < 1e-12);
        assert!(w[5] - 1.0 / 36.0 < 1e-12);
    }

    #[test]
    fn test_get_w_d3q19() {
        let vel_set_params = Parameters::test_default(3);

        let w = vel_set_params.get_w();

        assert!(w[0] - 1.0 / 3.0 < 1e-12);
        assert!(w[1] - 1.0 / 18.0 < 1e-12);
        assert!(w[7] - 1.0 / 36.0 < 1e-12);
    }

    #[test]
    fn test_get_q_bar_d2q9() {
        let vel_set_params = Parameters::test_default(2);

        let q_bar = vel_set_params.get_q_bar();

        assert_eq!(q_bar[0], 0);
        assert_eq!(q_bar[1], 3);
        assert_eq!(q_bar[5], 7);
    }

    #[test]
    fn test_get_q_bar_d3q19() {
        let vel_set_params = Parameters::test_default(3);

        let q_bar = vel_set_params.get_q_bar();

        assert_eq!(q_bar[0], 0);
        assert_eq!(q_bar[1], 2);
        assert_eq!(q_bar[7], 8);
        assert_eq!(q_bar[9], 10);
    }
}
