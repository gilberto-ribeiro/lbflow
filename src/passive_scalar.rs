pub mod bc;
mod io;
mod lattice;
mod node;
mod post;

use crate::cli;
use crate::prelude::*;
use bc::BoundaryCondition;
pub use lattice::Lattice;
pub use node::Node;

pub struct Parameters {
    pub scalar_name: String,
    pub tau_g: Float,
    pub initial_concentration: Vec<Float>,
    pub velocity_set: VelocitySet,
    pub boundary_conditions: Vec<(BoundaryFace, BoundaryCondition)>,
}

#[derive(Debug, Clone)]
pub struct Residuals {
    pub concentration: Float,
}

impl Residuals {
    pub fn new(concentration: Float) -> Self {
        Residuals { concentration }
    }
}

impl Residuals {
    pub fn get_concentration(&self) -> Float {
        self.concentration
    }

    pub fn set_concentration(&mut self, concentration: Float) {
        self.concentration = concentration;
    }
}

#[derive(Debug)]
pub struct ConversionFactor {
    pub tau_g: Float,
    pub diffusion_coefficient_conversion_factor: Float,
    pub diffusion_coefficient: Float,
    pub physical_diffusion_coefficient: Float,
}

impl ConversionFactor {
    pub fn new(tau_g: Float, delta_x: Float, delta_t: Float) -> Self {
        let diffusion_coefficient_conversion_factor = delta_x * delta_x / delta_t;
        let diffusion_coefficient = CS_2 * (tau_g - 0.5 * DELTA_T);
        let physical_diffusion_coefficient =
            diffusion_coefficient * diffusion_coefficient_conversion_factor;
        ConversionFactor {
            tau_g,
            diffusion_coefficient_conversion_factor,
            diffusion_coefficient,
            physical_diffusion_coefficient,
        }
    }

    pub fn from(
        momentum_params: &momentum::Parameters,
        passive_scalar_params: &Parameters,
    ) -> Self {
        let tau_g = passive_scalar_params.tau_g;
        let delta_x = momentum_params.delta_x;
        let delta_t = momentum_params.delta_t;
        ConversionFactor::new(tau_g, delta_x, delta_t)
    }
}

// ----------------------------------------------------------------------------- FUNCTIONS

pub fn run(
    config: Config,
    momentum_params: momentum::Parameters,
    passive_scalar_params: Parameters,
) {
    momentum::io::case_setup(&config, &momentum_params);

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
        cli::Mode::Post => post::vtk::post_vtk(config, momentum_parameters, passive_scalar_parameters),
        // cli::Mode::Post => println!("Post-processing is not implemented for passive scalar yet."),
    }
}
