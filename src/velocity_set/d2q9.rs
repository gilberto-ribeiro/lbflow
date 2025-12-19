use crate::BoundaryFace;
use crate::constants::Float;

pub(super) const D: usize = 2;

pub(super) const Q: usize = 9;

pub(super) const C: [[i32; D]; Q] = [
    [0, 0],
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1],
    [1, 1],
    [-1, 1],
    [-1, -1],
    [1, -1],
];

pub(super) const W: [Float; Q] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

pub(super) const Q_BAR: [usize; Q] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

const Q_WEST: [usize; 3] = [3, 6, 7];

const Q_EAST: [usize; 3] = [1, 5, 8];

const Q_SOUTH: [usize; 3] = [4, 7, 8];

const Q_NORTH: [usize; 3] = [2, 5, 6];

pub(super) const Q_FACES: [[usize; 3]; 4] = [Q_WEST, Q_EAST, Q_SOUTH, Q_NORTH];

pub(super) const MRT_MATRIX: [[Float; Q]; Q] = [
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [-4.0, -1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0],
    [4.0, -2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0],
    [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0],
    [0.0, -2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0],
    [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0],
    [0.0, 0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0],
    [0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0],
];

pub(super) const MRT_INVERSE_MATRIX: [[Float; Q]; Q] = [
    [
        1.0 / 9.0,
        -1.0 / 9.0,
        1.0 / 9.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        1.0 / 9.0,
        -1.0 / 36.0,
        -1.0 / 18.0,
        1.0 / 6.0,
        -1.0 / 6.0,
        0.0,
        0.0,
        0.25,
        0.0,
    ],
    [
        1.0 / 9.0,
        -1.0 / 36.0,
        -1.0 / 18.0,
        0.0,
        0.0,
        1.0 / 6.0,
        -1.0 / 6.0,
        -0.25,
        0.0,
    ],
    [
        1.0 / 9.0,
        -1.0 / 36.0,
        -1.0 / 18.0,
        -1.0 / 6.0,
        1.0 / 6.0,
        0.0,
        0.0,
        0.25,
        0.0,
    ],
    [
        1.0 / 9.0,
        -1.0 / 36.0,
        -1.0 / 18.0,
        0.0,
        0.0,
        -1.0 / 6.0,
        1.0 / 6.0,
        -0.25,
        0.0,
    ],
    [
        1.0 / 9.0,
        1.0 / 18.0,
        1.0 / 36.0,
        1.0 / 6.0,
        1.0 / 12.0,
        1.0 / 6.0,
        1.0 / 12.0,
        0.0,
        0.25,
    ],
    [
        1.0 / 9.0,
        1.0 / 18.0,
        1.0 / 36.0,
        -1.0 / 6.0,
        -1.0 / 12.0,
        1.0 / 6.0,
        1.0 / 12.0,
        0.0,
        -0.25,
    ],
    [
        1.0 / 9.0,
        1.0 / 18.0,
        1.0 / 36.0,
        -1.0 / 6.0,
        -1.0 / 12.0,
        -1.0 / 6.0,
        -1.0 / 12.0,
        0.0,
        0.25,
    ],
    [
        1.0 / 9.0,
        1.0 / 18.0,
        1.0 / 36.0,
        1.0 / 6.0,
        1.0 / 12.0,
        -1.0 / 6.0,
        -1.0 / 12.0,
        0.0,
        -0.25,
    ],
];

pub(super) fn velocity_computation(density: Float, f: Vec<Float>) -> Vec<Float> {
    vec![
        (1.0 / density) * (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]),
        (1.0 / density) * (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]),
    ]
}

pub(super) fn zou_he_bc_computation(
    boundary_face: &BoundaryFace,
    f: &[Float],
    density: &Option<Float>,
    velocity: &[Option<Float>],
) -> Vec<Float> {
    let velocity = &velocity[..D];
    let mut f = f.to_owned();
    match boundary_face {
        BoundaryFace::West => match (density, velocity) {
            (None, [Some(ux), Some(uy)]) => {
                let rho = 1.0 / (1.0 - ux) * (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7]));
                f[1] = f[3] + (2.0 / 3.0) * rho * ux;
                f[5] = f[7] - 0.5 * (f[2] - f[4]) + 0.5 * rho * uy + (1.0 / 6.0) * rho * ux;
                f[8] = f[6] + 0.5 * (f[2] - f[4]) - 0.5 * rho * uy + (1.0 / 6.0) * rho * ux;
            }
            (Some(rho), [None, Some(uy)]) => {
                let ux = 1.0 - 1.0 / rho * (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7]));
                f[1] = f[3] + (2.0 / 3.0) * rho * ux;
                f[5] = f[7] - 0.5 * (f[2] - f[4]) + 0.5 * rho * uy + (1.0 / 6.0) * rho * ux;
                f[8] = f[6] + 0.5 * (f[2] - f[4]) - 0.5 * rho * uy + (1.0 / 6.0) * rho * ux;
            }
            _ => {
                panic!("{}", crate::momentum::bc::ZOUHE_ERROR_MESSAGE);
            }
        },
        BoundaryFace::East => match (density, velocity) {
            (None, [Some(ux), Some(uy)]) => {
                let rho = 1.0 / (1.0 + ux) * (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8]));
                f[3] = f[1] - (2.0 / 3.0) * rho * ux;
                f[7] = f[5] + 0.5 * (f[2] - f[4]) - 0.5 * rho * uy - (1.0 / 6.0) * rho * ux;
                f[6] = f[8] - 0.5 * (f[2] - f[4]) + 0.5 * rho * uy - (1.0 / 6.0) * rho * ux;
            }
            (Some(rho), [None, Some(uy)]) => {
                let ux = 1.0 / rho * (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) - 1.0;
                f[3] = f[1] - (2.0 / 3.0) * rho * ux;
                f[7] = f[5] + 0.5 * (f[2] - f[4]) - 0.5 * rho * uy - (1.0 / 6.0) * rho * ux;
                f[6] = f[8] - 0.5 * (f[2] - f[4]) + 0.5 * rho * uy - (1.0 / 6.0) * rho * ux;
            }
            _ => {
                panic!("{}", crate::momentum::bc::ZOUHE_ERROR_MESSAGE);
            }
        },
        BoundaryFace::South => match (density, velocity) {
            (None, [Some(ux), Some(uy)]) => {
                let rho = 1.0 / (1.0 - uy) * (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8]));
                f[2] = f[4] + (2.0 / 3.0) * rho * uy;
                f[5] = f[7] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux + (1.0 / 6.0) * rho * uy;
                f[6] = f[8] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux + (1.0 / 6.0) * rho * uy;
            }
            (Some(rho), [Some(ux), None]) => {
                let uy = 1.0 - 1.0 / rho * (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8]));
                f[2] = f[4] + (2.0 / 3.0) * rho * uy;
                f[5] = f[7] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux + (1.0 / 6.0) * rho * uy;
                f[6] = f[8] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux + (1.0 / 6.0) * rho * uy;
            }
            _ => {
                panic!("{}", crate::momentum::bc::ZOUHE_ERROR_MESSAGE);
            }
        },
        BoundaryFace::North => match (density, velocity) {
            (None, [Some(ux), Some(uy)]) => {
                let rho = 1.0 / (1.0 + uy) * (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6]));
                f[4] = f[2] - (2.0 / 3.0) * rho * uy;
                f[7] = f[5] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux - (1.0 / 6.0) * rho * uy;
                f[8] = f[6] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux - (1.0 / 6.0) * rho * uy;
            }
            (Some(rho), [Some(ux), None]) => {
                let uy = 1.0 / rho * (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6])) - 1.0;
                f[4] = f[2] - (2.0 / 3.0) * rho * uy;
                f[7] = f[5] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux - (1.0 / 6.0) * rho * uy;
                f[8] = f[6] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux - (1.0 / 6.0) * rho * uy;
            }
            _ => {
                panic!("{}", crate::momentum::bc::ZOUHE_ERROR_MESSAGE);
            }
        },
        _ => panic!("{}", crate::momentum::bc::ZOUHE_ERROR_MESSAGE),
    }
    f
}

pub(super) fn mrt_equilibrium_moments_computation(
    density: Float,
    velocity: Vec<Float>,
) -> Vec<Float> {
    vec![
        density,
        density - 3.0 * density * (velocity[0] * velocity[0] + velocity[1] * velocity[1]),
        9.0 * density * velocity[0] * velocity[0] * velocity[1] * velocity[1]
            - 3.0 * density * (velocity[0] * velocity[0] + velocity[1] * velocity[1])
            + density,
        density * velocity[0],
        3.0 * density * velocity[0] * velocity[0] * velocity[0] - density * velocity[0],
        density * velocity[1],
        3.0 * density * velocity[1] * velocity[1] * velocity[1] - density * velocity[1],
        density * (velocity[0] * velocity[0] + velocity[1] * velocity[1]),
        density * velocity[0] * velocity[1],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_velocity_computation_d2q9() {
        let density = 1.0;
        let f = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];

        let velocity = velocity_computation(density, f.clone());

        let actual = velocity;
        let target = [-0.2, -0.6];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }

        let density = 0.5;

        let velocity = velocity_computation(density, f.clone());

        let actual = velocity;
        let target = [-0.4, -1.2];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }
}
