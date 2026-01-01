// ------------------------------------------------------------------------------- MODULES

pub mod bc;
pub(crate) mod io;
mod lattice;
pub(crate) mod multiphase;
mod node;
pub(crate) mod post;

// ------------------------------------------------------------------------------- IMPORTS

use crate::cli;
use crate::prelude_crate::*;
use bc::BoundaryCondition;
use clap::Parser;
pub(crate) use lattice::Lattice;
pub(crate) use node::Node;
pub(crate) use post::PostFunction;

// -------------------------------------------------------------------- STRUCT: Parameters

pub struct Parameters<'a> {
    pub n: Vec<usize>,
    pub collision_operator: CollisionOperator,
    pub force: Option<Box<dyn Fn(&Node) -> Vec<Float> + Send + Sync>>,
    pub multiphase_parameters: Option<multiphase::Parameters>,
    pub delta_x: Float,
    pub delta_t: Float,
    pub physical_density: Float,
    pub reference_pressure: Float,
    pub initial_density: InitialDensity<'a>,
    pub initial_velocity: InitialVelocity<'a>,
    pub velocity_set: VelocitySet,
    pub boundary_conditions: Vec<(BoundaryFace, BoundaryCondition)>,
    pub node_types: NodeTypes,
    pub post_functions: Option<Vec<PostFunction>>,
}

impl<'a> Default for Parameters<'a> {
    fn default() -> Self {
        Parameters {
            n: vec![10, 10],
            collision_operator: BGK(0.5),
            force: None,
            multiphase_parameters: None,
            delta_x: 0.01,
            delta_t: 0.01,
            physical_density: 998.0,
            reference_pressure: 101325.0,
            initial_density: InitialDensity::Uniform(1.0),
            initial_velocity: InitialVelocity::Uniform(vec![0.0, 0.0]),
            velocity_set: VelocitySet::D2Q9,
            node_types: NodeTypes::OnlyFluidNodes,
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

impl<'a> Parameters<'a> {
    fn _test_default(dim: usize) -> Self {
        match dim {
            2 => Default::default(),
            3 => Parameters {
                n: vec![10, 10, 10],
                initial_density: InitialDensity::Uniform(1.0),
                initial_velocity: InitialVelocity::Uniform(vec![0.0, 0.0, 0.0]),
                velocity_set: VelocitySet::D3Q19,
                node_types: NodeTypes::OnlyFluidNodes,
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
pub(super) struct Residuals {
    density: Float,
    velocity: Vec<Float>,
}

impl Residuals {
    pub fn new(density: Float, velocity: Vec<Float>) -> Self {
        Residuals { density, velocity }
    }
}

impl Residuals {
    pub(super) fn get_density(&self) -> Float {
        self.density
    }

    fn _set_density(&mut self, density: Float) {
        self.density = density;
    }

    pub(super) fn get_velocity(&self) -> &Vec<Float> {
        &self.velocity
    }

    fn _set_velocity(&mut self, velocity: Vec<Float>) {
        self.velocity = velocity;
    }
}

#[derive(Debug)]
pub(crate) struct ConversionFactor {
    pub(crate) collision_operator: Arc<CollisionOperator>,
    pub(crate) velocity_set: VelocitySet,
    pub(crate) delta_x: Float,
    pub(crate) delta_t: Float,
    pub(crate) length_conversion_factor: Float,
    pub(crate) time_conversion_factor: Float,
    pub(crate) density_conversion_factor: Float,
    pub(crate) velocity_conversion_factor: Float,
    pub(crate) pressure_conversion_factor: Float,
    pub(crate) viscosity_conversion_factor: Float,
    pub(crate) physical_density: Float,
    pub(crate) reference_pressure: Float,
    pub(crate) shear_viscosity: Float,
    pub(crate) bulk_viscosity: Float,
    pub(crate) physical_shear_viscosity: Float,
    pub(crate) physical_bulk_viscosity: Float,
}

impl ConversionFactor {
    fn new(
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
            MRT(relaxation_vector) => match velocity_set {
                D2Q5 => todo!("Implementing conversion factor for D2Q5."),
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
            },
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

    pub(crate) fn from(params: Parameters) -> Self {
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
        ConversionFactor::new(
            Arc::new(CollisionOperator::default()),
            VelocitySet::D2Q9,
            0.01,
            0.01,
            998.0,
            101325.0,
        )
    }
}

// ----------------------------------------------------------------------------- FUNCTIONS

fn run(cli_args: Cli, momentum_params: Parameters) {
    io::case_setup(&momentum_params);

    let lat = Lattice::new(cli_args, momentum_params);

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
    let cli_args = Cli::parse();

    let number_of_threads = cli_args.get_number_of_threads();
    crate::cli::init_global_pool(number_of_threads, cli_args.core_affinity);

    match &cli_args.command {
        cli::Command::Run { .. } => run(cli_args, momentum_params),
        cli::Command::Post { command } => match command {
            cli::PostCommand::Vtk { .. } => post::vtk::post_vtk(cli_args, momentum_params),
            cli::PostCommand::Unify { .. } => post::vtk::post_unify(cli_args, momentum_params),
        },
    }
}
