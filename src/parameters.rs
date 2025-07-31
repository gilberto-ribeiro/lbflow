pub mod functions;

use crate::BoundaryFace::*;
use crate::io::WriteDataMode;
use crate::momentum;
use crate::momentum::bc::BoundaryCondition::*;
use crate::parameters::functions::{fluid_and_solid_nodes_by_bounce_back_map, only_fluid_nodes};
use crate::simulation;
use crate::velocity_set::VelocitySet::*;

pub fn get_momentum_parameters() -> momentum::Parameters {
    momentum::Parameters {
        n: vec![200, 20, 20],
        tau: 0.9,
        physical_delta_x: 1.0,
        physical_delta_t: 1.0,
        velocity_set: D3Q15,
        node_types: only_fluid_nodes(vec![200, 20, 20]),
        boundary_conditions: vec![
            (West, AntiBounceBack { density: 1.05 }),
            (East, AntiBounceBack { density: 1.00 }),
            (South, NoSlip),
            (North, NoSlip),
            (Bottom, NoSlip),
            (Top, NoSlip),
        ],
    }
}

pub fn get_simulation_parameters() -> simulation::Parameters {
    simulation::Parameters {
        case_name: String::from("test_case"),
        data_write_mode: WriteDataMode::Frequency(100),
    }
}
