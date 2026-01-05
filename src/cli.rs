use clap::{Parser, Subcommand};
use core_affinity::{get_core_ids, set_for_current};
use std::num::NonZeroUsize;

pub(crate) type LbResult<T> = Result<T, Box<dyn std::error::Error>>;

const WRITE_DATA: usize = 100;
const MAX_ITERATIONS: usize = 100_000;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None, propagate_version = true, subcommand_required = true, arg_required_else_help = true)]
pub(crate) struct Cli {
    /// Number of Rayon threads (min = 1)
    #[arg(
        short,
        long = "num-threads",
        alias = "threads",
        value_name = "NTHREADS",
        default_value = "1",
        global = true
    )]
    pub(crate) num_threads: NonZeroUsize,
    /// Enable core pinning (affinity)
    #[arg(short, long, alias = "affinity", global = true)]
    pub(crate) pin: bool,
    /// Start index (within eligible cores list) for pinning
    #[arg(
        long = "core-start",
        value_name = "IDX",
        default_value_t = 0,
        global = true
    )]
    pub(crate) core_start: usize,
    /// Number of cores to use starting from core-start (optional)
    #[arg(long = "core-count", value_name = "N", global = true)]
    pub(crate) core_count: Option<usize>,
    /// Explicit CPU list, e.g.: "0-3,8,10-12" (overrides core-start/core-count)
    #[arg(long = "cores", value_name = "LIST", global = true)]
    pub(crate) cores_list: Option<String>,
    /// Freeze the momentum field calculation and post-processing
    #[arg(short, long = "freeze-momentum", alias = "freeze", global = true)]
    freeze_momentum: bool,
    #[command(subcommand)]
    pub(crate) command: Command,
}

impl Default for Cli {
    fn default() -> Self {
        Self {
            num_threads: NonZeroUsize::new(1).unwrap(),
            pin: false,
            core_start: 0,
            core_count: None,
            cores_list: None,
            freeze_momentum: false,
            command: Command::Run {
                write_data: WRITE_DATA,
                max_iterations: MAX_ITERATIONS,
            },
        }
    }
}

impl Cli {
    pub(crate) fn get_num_threads(&self) -> usize {
        usize::from(self.num_threads)
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
        self.freeze_momentum
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

pub(crate) fn init_global_pool(
    num_threads: usize,
    pin: bool,
    core_start: usize,
    core_count: Option<usize>,
    cores_list: Option<&str>,
) {
    if pin {
        let eligible = get_core_ids().expect("Failed to query system CPU cores");

        // Dynamic default:
        // If --pin and the user did NOT provide --core-count and did NOT provide --cores,
        // then core_count defaults to num_threads (1:1 thread-to-core).
        let core_count_eff = if cores_list.is_none() && core_count.is_none() {
            Some(num_threads)
        } else {
            core_count
        };

        let selected = select_cores(eligible, core_start, core_count_eff, cores_list);

        if num_threads > selected.len() {
            eprintln!(
                "Warning: --threads ({}) is greater than the number of selected cores ({}). \
Some threads will be pinned to the same core (oversubscription). \
Consider increasing --core-count/--cores or reducing --threads.",
                num_threads,
                selected.len()
            );
        }

        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .start_handler(move |idx| {
                let core = selected[idx % selected.len()];
                let _ = set_for_current(core);
            })
            .build_global()
            .expect("Rayon global thread pool was already initialized");
    } else {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .expect("Rayon global thread pool was already initialized");
    };
}

fn select_cores(
    mut eligible: Vec<core_affinity::CoreId>,
    core_start: usize,
    core_count: Option<usize>,
    explicit_cpus: Option<&str>,
) -> Vec<core_affinity::CoreId> {
    if let Some(list_str) = explicit_cpus {
        let cpus = parse_cpu_list(list_str).expect("Failed to parse --cores");
        let set: std::collections::HashSet<usize> = cpus.into_iter().collect();
        eligible.retain(|c| set.contains(&c.id));

        if eligible.is_empty() {
            panic!("--cores resulted in an empty eligible CPU set (check your list).");
        }
        return eligible;
    }

    let start = core_start.min(eligible.len());
    let end = match core_count {
        Some(n) => (start + n).min(eligible.len()),
        None => eligible.len(),
    };

    let selected = if start < end {
        eligible[start..end].to_vec()
    } else {
        eligible
    };

    if selected.is_empty() {
        panic!("No eligible CPU cores available for pinning.");
    }

    selected
}

fn parse_cpu_list(s: &str) -> Result<Vec<usize>, String> {
    let mut out = Vec::new();

    for part in s.split(',').map(str::trim).filter(|p| !p.is_empty()) {
        if let Some((a, b)) = part.split_once('-') {
            let start: usize = a
                .trim()
                .parse()
                .map_err(|_| format!("Invalid CPU id: {a}"))?;
            let end: usize = b
                .trim()
                .parse()
                .map_err(|_| format!("Invalid CPU id: {b}"))?;
            if end < start {
                return Err(format!("Invalid CPU range: {part}"));
            }
            out.extend(start..=end);
        } else {
            let cpu: usize = part
                .parse()
                .map_err(|_| format!("Invalid CPU id: {part}"))?;
            out.push(cpu);
        }
    }

    out.sort_unstable();
    out.dedup();
    Ok(out)
}
