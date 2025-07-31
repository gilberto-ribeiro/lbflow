use crate::constants::Float;

pub const D: usize = 3;

pub const Q: usize = 27;

pub const C: [[i32; D]; Q] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, -1],
    [1, -1, 0],
    [-1, 1, 0],
    [1, 0, -1],
    [-1, 0, 1],
    [0, 1, -1],
    [0, -1, 1],
    [1, 1, 1],
    [-1, -1, -1],
    [1, 1, -1],
    [-1, -1, 1],
    [1, -1, 1],
    [-1, 1, -1],
    [-1, 1, 1],
    [1, -1, -1],
];

pub const W: [Float; Q] = [
    8.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
];

pub const Q_BAR: [usize; Q] = [
    0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26,
    25,
];

const Q_WEST: [usize; 9] = [2, 8, 10, 14, 16, 20, 22, 24, 25];

const Q_EAST: [usize; 9] = [1, 7, 9, 13, 15, 19, 21, 23, 26];

const Q_SOUTH: [usize; 9] = [4, 8, 12, 13, 18, 20, 22, 23, 26];

const Q_NORTH: [usize; 9] = [3, 7, 11, 14, 17, 19, 21, 24, 25];

const Q_BOTTOM: [usize; 9] = [6, 10, 12, 15, 17, 20, 21, 24, 26];

const Q_TOP: [usize; 9] = [5, 9, 11, 16, 18, 19, 22, 23, 25];

pub const Q_FACES: [[usize; 9]; 6] = [Q_WEST, Q_EAST, Q_SOUTH, Q_NORTH, Q_BOTTOM, Q_TOP];

/// # Examples:
/// ```
/// # use lbflow::velocity_set::d3q27;
/// let density = 1.0;
/// //             0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26
/// let f = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
/// let velocity = d3q27::velocity_computation(density, f.clone());
/// let actual = velocity;
/// let target = vec![-0.1, 0.1, 0.3];
/// for (a, b) in actual.iter().zip(target.iter()) {
///     assert!((a - b).abs() < 1e-12);
/// }
///
/// let density = 0.5;
/// let velocity = d3q27::velocity_computation(density, f.clone());
/// let actual = velocity;
/// let target = vec![-0.2, 0.2, 0.6];
/// for (a, b) in actual.iter().zip(target.iter()) {
///     assert!((a - b).abs() < 1e-12);
/// }
/// ```
pub fn velocity_computation(density: Float, f: Vec<Float>) -> Vec<Float> {
    vec![
        (1.0 / density)
            * (f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[13] - f[14] + f[15] - f[16] + f[19]
                - f[20]
                + f[21]
                - f[22]
                + f[23]
                - f[24]
                - f[25]
                + f[26]),
        (1.0 / density)
            * (f[3] - f[4] + f[7] - f[8] + f[11] - f[12] - f[13] + f[14] + f[17] - f[18] + f[19]
                - f[20]
                + f[21]
                - f[22]
                - f[23]
                + f[24]
                + f[25]
                - f[26]),
        (1.0 / density)
            * (f[5] - f[6] + f[9] - f[10] + f[11] - f[12] - f[15] + f[16] - f[17] + f[18] + f[19]
                - f[20]
                - f[21]
                + f[22]
                + f[23]
                - f[24]
                + f[25]
                - f[26]),
    ]
}
