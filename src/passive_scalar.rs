pub(crate) mod adsorption;
pub mod bc;
pub(crate) mod io;
mod lattice;
mod node;
mod post;

use crate::cli;
use crate::prelude_crate::*;
use bc::BoundaryCondition;
use bc::InnerBoundaryCondition;
pub(crate) use lattice::Lattice;
pub(crate) use node::Node;

pub struct Parameters<'a> {
    pub scalar_name: &'a str,
    pub collision_operator: CollisionOperator,
    pub source_value: Option<Box<dyn Fn(&Node) -> Float + Send + Sync>>,
    pub initial_scalar_value: InitialScalarValue<'a>,
    pub velocity_set: VelocitySet,
    pub boundary_conditions: Vec<(BoundaryFace, BoundaryCondition)>,
    pub inner_boundary_condition: InnerBoundaryCondition,
    pub adsorption_parameters: Option<adsorption::Parameters>,
}

#[derive(Debug, Clone)]
struct Residuals {
    scalar_value: Float,
}

impl Residuals {
    fn new(scalar_value: Float) -> Self {
        Residuals { scalar_value }
    }
}

impl Residuals {
    fn get_scalar_value(&self) -> Float {
        self.scalar_value
    }

    fn _set_scalar_value(&mut self, scalar_value: Float) {
        self.scalar_value = scalar_value;
    }
}

#[derive(Debug)]
pub(crate) struct _ConversionFactor {
    pub(crate) tau_g: Float,
    pub(crate) diffusion_coefficient_conversion_factor: Float,
    pub(crate) diffusion_coefficient: Float,
    pub(crate) physical_diffusion_coefficient: Float,
}

impl _ConversionFactor {
    fn _new(tau_g: Float, delta_x: Float, delta_t: Float) -> Self {
        let diffusion_coefficient_conversion_factor = delta_x * delta_x / delta_t;
        let diffusion_coefficient = CS_2 * (tau_g - 0.5 * DELTA_T);
        let physical_diffusion_coefficient =
            diffusion_coefficient * diffusion_coefficient_conversion_factor;
        _ConversionFactor {
            tau_g,
            diffusion_coefficient_conversion_factor,
            diffusion_coefficient,
            physical_diffusion_coefficient,
        }
    }
}

// ----------------------------------------------------------------------------- FUNCTIONS

fn run(config: Config, momentum_params: momentum::Parameters, passive_scalar_params: Parameters) {
    momentum::io::case_setup(&momentum_params);

    let m_lat = Arc::new(momentum::Lattice::new(config, momentum_params));
    let ps_lat = Lattice::new(passive_scalar_params, Arc::clone(&m_lat));

    m_lat.write_coordinates().unwrap_or_else(|e| {
        eprintln! {"Error while writing the coordinates file: {e}"};
        std::process::exit(1);
    });

    m_lat.initialize_nodes();
    ps_lat.initialize_nodes();
    loop {
        m_lat.main_steps();
        m_lat.compute_lattice_residuals();

        ps_lat.main_steps();
        ps_lat.compute_lattice_residuals();

        m_lat.write_data();
        ps_lat.write_data();

        m_lat.compute_post_processing();

        crate::io::print_residuals(&ps_lat.get_residuals_info());
        crate::io::write_residuals(&ps_lat.get_residuals_info()).unwrap_or_else(|e| {
            eprintln! {"Error while writing the residuals file: {e}"};
            std::process::exit(1);
        });

        if ps_lat.stop_condition() {
            break;
        };

        m_lat.update_shallow_nodes();
        ps_lat.update_shallow_nodes();

        m_lat.next_time_step();
    }
}

pub fn solve(momentum_parameters: momentum::Parameters, passive_scalar_parameters: Parameters) {
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
        cli::Mode::Run => run(config, momentum_parameters, passive_scalar_parameters),
        cli::Mode::PostVTK => {
            post::vtk::post_vtk(config, momentum_parameters, passive_scalar_parameters)
        }
        cli::Mode::PostUnify => {
            post::vtk::post_unify(config, momentum_parameters, passive_scalar_parameters)
        }
    }
}

fn run_vec(
    config: Config,
    momentum_params: momentum::Parameters,
    passive_scalar_params_vec: Vec<Parameters>,
) {
    momentum::io::case_setup(&momentum_params);

    let m_lat = Arc::new(momentum::Lattice::new(config, momentum_params));
    let ps_lat_vec = lattice::LatticeVec::new(passive_scalar_params_vec, Arc::clone(&m_lat));

    m_lat.write_coordinates().unwrap_or_else(|e| {
        eprintln! {"Error while writing the coordinates file: {e}"};
        std::process::exit(1);
    });

    m_lat.initialize_nodes();
    ps_lat_vec.initialize_nodes();
    loop {
        m_lat.main_steps();
        m_lat.compute_lattice_residuals();

        ps_lat_vec.main_steps();
        ps_lat_vec.compute_lattice_residuals();

        m_lat.write_data();
        ps_lat_vec.write_data();

        m_lat.compute_post_processing();

        crate::io::print_residuals(&ps_lat_vec.get_residuals_info());
        crate::io::write_residuals(&ps_lat_vec.get_residuals_info()).unwrap_or_else(|e| {
            eprintln! {"Error while writing the residuals file: {e}"};
            std::process::exit(1);
        });

        if ps_lat_vec.stop_condition() {
            break;
        };

        m_lat.update_shallow_nodes();
        ps_lat_vec.update_shallow_nodes();

        m_lat.next_time_step();
    }
}

pub fn solve_vec(
    momentum_parameters: momentum::Parameters,
    passive_scalar_parameters_vec: Vec<Parameters>,
) {
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
        cli::Mode::Run => run_vec(config, momentum_parameters, passive_scalar_parameters_vec),
        cli::Mode::PostVTK => {
            post::vtk::post_vtk_vec(config, momentum_parameters, passive_scalar_parameters_vec)
        }
        cli::Mode::PostUnify => {
            post::vtk::post_unify_vec(config, momentum_parameters, passive_scalar_parameters_vec)
        }
    }
}
