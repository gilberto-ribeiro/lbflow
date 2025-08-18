// ------------------------------------------------------------------------------- MODULES

pub mod bc;
mod io;
mod lattice;
mod node;
pub mod post;

// ------------------------------------------------------------------------------- IMPORTS

use crate::cli;
use crate::prelude::*;
use bc::BoundaryCondition;
pub use lattice::Lattice;
pub use node::{Node, ShallowNode};
pub use post::{PostFunction, PostResult};

// -------------------------------------------------------------------- STRUCT: Parameters

pub struct Parameters {
    pub n: Vec<usize>,
    pub tau: Float,
    pub delta_x: Float,
    pub delta_t: Float,
    pub physical_density: Float,
    pub reference_pressure: Float,
    pub initial_density: Vec<Float>,
    pub initial_velocity: Vec<Vec<Float>>,
    pub velocity_set: VelocitySet,
    pub boundary_conditions: Vec<(BoundaryFace, BoundaryCondition)>,
    pub node_types: Vec<NodeType>,
    pub post_functions: Option<Vec<PostFunction>>,
}

impl Default for Parameters {
    fn default() -> Self {
        Parameters {
            n: vec![10, 10],
            tau: 0.5,
            delta_x: 0.01,
            delta_t: 0.01,
            physical_density: 998.0,
            reference_pressure: 101325.0,
            initial_density: functions::uniform_density(1.0, vec![10, 10]),
            initial_velocity: functions::uniform_velocity(vec![0.0, 0.0], vec![10, 10]),
            velocity_set: VelocitySet::D2Q9,
            node_types: functions::only_fluid_nodes(vec![10, 10]),
            boundary_conditions: vec![
                (BoundaryFace::West, BoundaryCondition::NoSlip),
                (BoundaryFace::East, BoundaryCondition::NoSlip),
                (BoundaryFace::South, BoundaryCondition::NoSlip),
                (
                    BoundaryFace::North,
                    BoundaryCondition::BounceBack {
                        density: 1.0,
                        velocity: vec![0.1, 0.0],
                    },
                ),
            ],
            post_functions: None,
        }
    }
}

impl Parameters {
    pub fn test_default(dim: usize) -> Self {
        match dim {
            2 => Default::default(),
            3 => Parameters {
                n: vec![10, 10, 10],
                initial_density: functions::uniform_density(1.0, vec![10, 10, 10]),
                initial_velocity: functions::uniform_velocity(
                    vec![0.0, 0.0, 0.0],
                    vec![10, 10, 10],
                ),
                velocity_set: VelocitySet::D3Q27,
                node_types: functions::only_fluid_nodes(vec![10, 10, 10]),
                boundary_conditions: vec![
                    (BoundaryFace::West, BoundaryCondition::NoSlip),
                    (BoundaryFace::East, BoundaryCondition::NoSlip),
                    (BoundaryFace::South, BoundaryCondition::NoSlip),
                    (
                        BoundaryFace::North,
                        BoundaryCondition::BounceBack {
                            density: 1.0,
                            velocity: vec![0.1, 0.0, 0.0],
                        },
                    ),
                    (BoundaryFace::Bottom, BoundaryCondition::NoSlip),
                    (BoundaryFace::Top, BoundaryCondition::NoSlip),
                ],
                ..Default::default()
            },
            _ => panic!("Unsupported dimension: {dim}"),
        }
    }
}

// --------------------------------------------------------------------- STRUCT: Residuals

#[derive(Debug, Clone)]
pub struct Residuals {
    pub density: Float,
    pub velocity: Vec<Float>,
}

impl Residuals {
    pub fn new(density: Float, velocity: Vec<Float>) -> Self {
        Residuals { density, velocity }
    }
}

impl Residuals {
    pub fn get_density(&self) -> Float {
        self.density
    }

    pub fn set_density(&mut self, density: Float) {
        self.density = density;
    }

    pub fn get_velocity(&self) -> &Vec<Float> {
        &self.velocity
    }

    pub fn set_velocity(&mut self, velocity: Vec<Float>) {
        self.velocity = velocity;
    }
}

#[derive(Debug)]
pub struct ConversionFactor {
    pub tau: Float,
    pub delta_x: Float,
    pub delta_t: Float,
    pub length_conversion_factor: Float,
    pub time_conversion_factor: Float,
    pub density_conversion_factor: Float,
    pub velocity_conversion_factor: Float,
    pub pressure_conversion_factor: Float,
    pub viscosity_conversion_factor: Float,
    pub physical_density: Float,
    pub reference_pressure: Float,
    pub viscosity: Float,
    pub physical_viscosity: Float,
}

impl ConversionFactor {
    pub fn new(
        tau: Float,
        delta_x: Float,
        delta_t: Float,
        physical_density: Float,
        reference_pressure: Float,
    ) -> Self {
        let length_conversion_factor = delta_x;
        let time_conversion_factor = delta_t;
        let density_conversion_factor = physical_density;
        let velocity_conversion_factor = delta_x / delta_t;
        let pressure_conversion_factor =
            density_conversion_factor * velocity_conversion_factor * velocity_conversion_factor;
        let viscosity_conversion_factor = delta_x * delta_x / delta_t;
        let viscosity = CS_2 * (tau - 0.5);
        let physical_viscosity = viscosity_conversion_factor * viscosity;
        ConversionFactor {
            tau,
            delta_x,
            delta_t,
            length_conversion_factor,
            time_conversion_factor,
            density_conversion_factor,
            velocity_conversion_factor,
            pressure_conversion_factor,
            viscosity_conversion_factor,
            physical_density,
            reference_pressure,
            viscosity,
            physical_viscosity,
        }
    }

    pub fn from(params: &Parameters) -> Self {
        ConversionFactor::new(
            params.tau,
            params.delta_x,
            params.delta_t,
            params.physical_density,
            params.reference_pressure,
        )
    }
}

impl Default for ConversionFactor {
    fn default() -> Self {
        ConversionFactor::new(0.5, 0.01, 0.01, 998.0, 101325.0)
    }
}

// ----------------------------------------------------------------------------- FUNCTIONS

pub fn run(config: Config, momentum_parameters: Parameters) {
    io::case_setup(&config, &momentum_parameters);
    let lattice = Lattice::new(config, momentum_parameters);
    lattice.write_coordinates().unwrap_or_else(|e| {
        eprintln! {"Error while writing the coordinates file: {e}"};
        std::process::exit(1);
    });
    lattice.initialize_nodes();
    loop {
        lattice.update_density_and_velocity_step();
        lattice.equilibrium_step();
        lattice.bgk_collision_step();
        lattice.streaming_step();
        lattice.inner_bounce_back_step();
        lattice.boundary_conditions_step();
        lattice.write_data();
        lattice.compute_post_processing();
        lattice.compute_lattice_residuals();
        lattice.print_residuals();
        lattice.write_residuals().unwrap_or_else(|e| {
            eprintln! {"Error while writing the residuals file: {e}"};
            std::process::exit(1);
        });
        if lattice.stop_condition() {
            break;
        };
        lattice.update_shallow_nodes();
        lattice.next_time_step();
    }
}

pub fn load(momentum_parameters: momentum::Parameters) {
    let config = match cli::get_args().and_then(|matches| cli::parse_matches(&matches)) {
        Ok(cfg) => cfg,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    };

    let number_of_threads = usize::from(config.number_of_threads);
    crate::cli::init_global_pool(number_of_threads, config.core_affinity);

    match config.mode {
        cli::Mode::Run => run(config, momentum_parameters),
        cli::Mode::Post => post::vtk::post_vtk(config, momentum_parameters),
    }
}
