use crate::constants::Float;

pub const D: usize = 3;

pub const Q: usize = 15;

pub const C: [[i32; D]; Q] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
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
    2.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
];

pub const Q_BAR: [usize; Q] = [0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13];

const Q_WEST: [usize; 5] = [2, 8, 10, 12, 13];

const Q_EAST: [usize; 5] = [1, 7, 9, 11, 14];

const Q_SOUTH: [usize; 5] = [4, 8, 10, 11, 14];

const Q_NORTH: [usize; 5] = [3, 7, 9, 12, 13];

const Q_BOTTOM: [usize; 5] = [6, 8, 9, 12, 14];

const Q_TOP: [usize; 5] = [5, 7, 10, 11, 13];

pub const Q_FACES: [[usize; 5]; 6] = [Q_WEST, Q_EAST, Q_SOUTH, Q_NORTH, Q_BOTTOM, Q_TOP];
