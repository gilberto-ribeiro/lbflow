use clap::{Parser, Subcommand};
use core_affinity::{get_core_ids, set_for_current};
use std::num::NonZeroUsize;

pub(crate) type LbResult<T> = Result<T, Box<dyn std::error::Error>>;

const WRITE_DATA: usize = 100;
const MAX_ITERATIONS: usize = 100_000;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None, propagate_version = true, subcommand_required = true, arg_required_else_help = true)]
pub(crate) struct Cli {
    /// The number of threads used (min = 1)
    #[arg(
        short,
        long = "num-threads",
        value_name = "NTHREADS",
        default_value = "1",
        global = true
    )]
    pub(crate) number_of_threads: NonZeroUsize,
    /// Set the core affinity
    #[arg(short = 'a', long = "affinity", global = true)]
    pub(crate) core_affinity: bool,
    #[command(subcommand)]
    pub(crate) command: Command,
}

impl Default for Cli {
    fn default() -> Self {
        Self {
            number_of_threads: NonZeroUsize::new(1).unwrap(),
            core_affinity: false,
            command: Command::Run {
                write_data: WRITE_DATA,
                max_iterations: MAX_ITERATIONS,
                freeze_momentum: false,
            },
        }
    }
}

impl Cli {
    pub(crate) fn get_number_of_threads(&self) -> usize {
        usize::from(self.number_of_threads)
    }

    pub(crate) fn get_write_data_frequency(&self) -> usize {
        match self.command {
            Command::Run { write_data, .. } => write_data,
            Command::Post { .. } => WRITE_DATA,
        }
    }

    pub(crate) fn get_max_iterations(&self) -> usize {
        match self.command {
            Command::Run { max_iterations, .. } => max_iterations,
            Command::Post { .. } => MAX_ITERATIONS,
        }
    }

    pub(crate) fn get_freeze_momentum(&self) -> bool {
        match self.command {
            Command::Run {
                freeze_momentum, ..
            } => freeze_momentum,
            Command::Post { .. } => false,
        }
    }
}

#[derive(Subcommand, Debug)]
pub(crate) enum Command {
    /// Run the simulation
    Run {
        /// The frequency wich data is written
        #[arg(
            short,
            long = "write-data",
            value_name = "FREQUENCY",
            default_value_t = WRITE_DATA,
        )]
        write_data: usize,
        /// The maximum number of iterations
        #[arg(
            long = "max-iterations",
            value_name = "ITER",
            default_value_t = MAX_ITERATIONS,
        )]
        max_iterations: usize,
        /// Freeze the momentum field calculation
        #[arg(short, long = "freeze")]
        freeze_momentum: bool,
    },
    /// Post-process the simulation data
    Post {
        #[command(subcommand)]
        command: PostCommand,
    },
}

impl Default for Command {
    fn default() -> Self {
        Self::Run {
            write_data: 100,
            max_iterations: 100_000,
            freeze_momentum: false,
        }
    }
}

#[derive(Subcommand, Debug)]
pub(crate) enum PostCommand {
    /// Unify the parallel output files into a single file for each time step
    Unify {
        /// Keep the intermediate files after unification
        #[arg(short, long)]
        keep: bool,
    },
    /// Write VTK files for visualization in Paraview
    Vtk {
        /// Write node_type.vtk file
        #[arg(long = "node-type")]
        node_type: bool,
        /// If the data is written in physical units
        #[arg(long = "physical")]
        physical_data: bool,
    },
}

impl Default for PostCommand {
    fn default() -> Self {
        Self::Vtk {
            node_type: false,
            physical_data: false,
        }
    }
}

pub(crate) fn init_global_pool(num_threads: usize, pin_all_cores: bool) {
    if pin_all_cores {
        let cores = get_core_ids().expect("listar cores do sistema");
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .start_handler(move |idx| {
                if pin_all_cores {
                    let core = cores[idx % cores.len()];
                    let _ = set_for_current(core);
                }
            })
            .build_global()
            .expect("pool global já foi criada antes?");
    } else {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .expect("pool global já foi criada antes?");
    };
}
