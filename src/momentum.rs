// ------------------------------------------------------------------------------- MODULES

pub mod bc;
pub mod io;
pub mod lattice;
mod node;
pub mod post;

// ------------------------------------------------------------------------------- IMPORTS

use crate::cli;
use crate::prelude::*;
use crate::velocity_set;
use bc::BoundaryCondition;
pub use lattice::Lattice;
pub use node::{Node, ShallowNode};
pub use post::{PostFunction, PostResult};

// -------------------------------------------------------------------- STRUCT: Parameters

pub struct Parameters {
    pub n: Vec<usize>,
    pub collision_operator: CollisionOperator,
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
            collision_operator: BGK(0.5),
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
    pub collision_operator: Arc<CollisionOperator>,
    pub velocity_set: VelocitySet,
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
    pub shear_viscosity: Float,
    pub bulk_viscosity: Float,
    pub physical_shear_viscosity: Float,
    pub physical_bulk_viscosity: Float,
}

impl ConversionFactor {
    pub fn new(
        collision_operator: Arc<CollisionOperator>,
        velocity_set: VelocitySet,
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
        let (shear_viscosity, bulk_viscosity) = match collision_operator.as_ref() {
            BGK(tau) => {
                let shear_viscosity = CS_2 * (tau - 0.5 * DELTA_T);
                let bulk_viscosity = 2.0 / 3.0 * shear_viscosity;
                (shear_viscosity, bulk_viscosity)
            }
            TRT(omega_plus, _) => {
                let shear_viscosity = CS_2 * (1.0 / (omega_plus * DELTA_T) - 0.5);
                let bulk_viscosity = 2.0 / 3.0 * shear_viscosity; // Confirmar
                (shear_viscosity, bulk_viscosity)
            }
            MRT(relaxation_vector) => {
                match velocity_set {
                    D2Q9 => {
                        let omega_nu = relaxation_vector[7];
                        let omega_e = relaxation_vector[1];
                        let shear_viscosity = CS_2 * (1.0 / omega_nu - 0.5);
                        let bulk_viscosity = CS_2 * (1.0 / omega_e - 0.5) - shear_viscosity / 3.0;
                        (shear_viscosity, bulk_viscosity)
                    }
                    D3Q15 => {
                        let omega_nu = relaxation_vector[9];
                        let omega_e = relaxation_vector[1];
                        let shear_viscosity = CS_2 * (1.0 / omega_nu - 0.5 * DELTA_T);
                        let bulk_viscosity = 2.0 / 3.0 * CS_2 * (1.0 / omega_e - 0.5 * DELTA_T);
                        (shear_viscosity, bulk_viscosity)
                    }
                    D3Q19 => {
                        let omega_nu = relaxation_vector[9];
                        let omega_e = relaxation_vector[1];
                        let shear_viscosity = CS_2 * (1.0 / omega_nu - 0.5 * DELTA_T);
                        let bulk_viscosity = 2.0 / 3.0 * CS_2 * (1.0 / omega_e - 0.5 * DELTA_T);
                        (shear_viscosity, bulk_viscosity)
                    }
                    D3Q27 => {
                        let omega_nu = relaxation_vector[9];
                        let omega_e = relaxation_vector[1];
                        let shear_viscosity = CS_2 * (1.0 / omega_nu - 0.5 * DELTA_T);
                        let bulk_viscosity = 2.0 / 3.0 * CS_2 * (1.0 / omega_e - 0.5 * DELTA_T);
                        (shear_viscosity, bulk_viscosity)
                    }
                }
            }
        };
        let physical_shear_viscosity = viscosity_conversion_factor * shear_viscosity;
        let physical_bulk_viscosity = viscosity_conversion_factor * bulk_viscosity;
        ConversionFactor {
            collision_operator,
            velocity_set,
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
            shear_viscosity,
            bulk_viscosity,
            physical_shear_viscosity,
            physical_bulk_viscosity,
        }
    }

    pub fn from(params: Parameters) -> Self {
        ConversionFactor::new(
            Arc::new(params.collision_operator),
            params.velocity_set,
            params.delta_x,
            params.delta_t,
            params.physical_density,
            params.reference_pressure,
        )
    }
}

impl Default for ConversionFactor {
    fn default() -> Self {
        ConversionFactor::new(Arc::new(CollisionOperator::default()), VelocitySet::D2Q9, 0.01, 0.01, 998.0, 101325.0)
    }
}

// ----------------------------------------------------------------------------- FUNCTIONS

pub fn run(config: Config, momentum_params: Parameters) {
    io::case_setup(&momentum_params);

    let lat = Lattice::new(config, momentum_params);

    lat.write_coordinates().unwrap_or_else(|e| {
        eprintln! {"Error while writing the coordinates file: {e}"};
        std::process::exit(1);
    });

    lat.initialize_nodes();
    loop {
        lat.main_steps();
        lat.compute_lattice_residuals();

        lat.write_data();
        lat.compute_post_processing();

        crate::io::print_residuals(&lat.get_residuals_info());
        crate::io::write_residuals(&lat.get_residuals_info()).unwrap_or_else(|e| {
            eprintln! {"Error while writing the residuals file: {e}"};
            std::process::exit(1);
        });

        if lat.stop_condition() {
            break;
        };

        lat.update_shallow_nodes();

        lat.next_time_step();
    }
}

pub fn solve(momentum_params: momentum::Parameters) {
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
        cli::Mode::Run => run(config, momentum_params),
        cli::Mode::Post => post::vtk::post_vtk(config, momentum_params),
    }
}
