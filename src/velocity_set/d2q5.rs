use crate::constants::Float;

pub(super) const D: usize = 2;

pub(super) const Q: usize = 5;

pub(super) const C: [[i32; D]; Q] = [[0, 0], [1, 0], [-1, 0], [0, 1], [0, -1]];

pub(super) const W: [Float; Q] = [2.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];

pub(super) const Q_BAR: [usize; Q] = [0, 2, 1, 4, 3];

const Q_WEST: [usize; 1] = [2];

const Q_EAST: [usize; 1] = [1];

const Q_SOUTH: [usize; 1] = [4];

const Q_NORTH: [usize; 1] = [3];

pub(super) const Q_FACES: [[usize; 1]; 4] = [Q_WEST, Q_EAST, Q_SOUTH, Q_NORTH];
