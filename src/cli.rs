use crate::prelude_crate::*;
use clap::{Arg, Command};
use core_affinity::{get_core_ids, set_for_current};
use std::num::{NonZero, NonZeroUsize};

pub(crate) type LbResult<T> = Result<T, Box<dyn std::error::Error>>;

#[derive(Debug)]
pub(crate) struct Config {
    pub(crate) mode: Mode,
    pub(crate) number_of_threads: NonZeroUsize,
    pub(crate) core_affinity: bool,
    pub(crate) write_data: Option<WriteDataMode>,
    pub(crate) max_iterations: Option<usize>,
    pub(crate) freeze_momentum: bool,
    pub(crate) node_type: bool,
    pub(crate) physical_data: bool,
    pub(crate) keep: bool,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            mode: Mode::Run,
            number_of_threads: NonZero::new(1).unwrap(),
            core_affinity: false,
            write_data: None,
            max_iterations: None,
            freeze_momentum: false,
            node_type: false,
            physical_data: false,
            keep: false,
        }
    }
}

impl Config {
    pub(crate) fn get_number_of_threads(&self) -> usize {
        usize::from(self.number_of_threads)
    }

    pub(crate) fn get_write_data_mode(&self) -> &WriteDataMode {
        self.write_data.as_ref().unwrap()
    }

    pub(crate) fn get_max_iterations(&self) -> usize {
        self.max_iterations.unwrap()
    }
}

#[derive(Debug)]
pub(crate) enum Mode {
    Run,
    PostVTK,
    PostUnify,
}

pub(crate) fn get_args() -> LbResult<clap::ArgMatches> {
    let matches = clap::command!()
        .propagate_version(true)
        .subcommand_required(true)
        .arg_required_else_help(true)
        .arg(
            Arg::new("number_of_threads")
                .short('n')
                .long("num-threads")
                .value_name("NTHREADS")
                .help("The number of threads used (min = 1)")
                .value_parser(clap::value_parser!(NonZeroUsize))
                .default_value("1")
                .global(true),
        )
        .arg(
            Arg::new("core_affinity")
                .long("affinity")
                .help("Set the core affinity")
                .action(clap::ArgAction::SetTrue)
                .global(true),
        )
        .subcommand(
            Command::new("run")
                .about("Run the simulation")
                .arg(
                    Arg::new("write_data")
                        .short('w')
                        .long("write-data")
                        .value_name("FREQUENCY")
                        .help("The frequency wich data is written")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("100"),
                )
                .arg(
                    Arg::new("max_iterations")
                        .long("max-iterations")
                        .value_name("ITER")
                        .help("The maximum number of iterations")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("100000"),
                )
                .arg(
                    Arg::new("freeze_momentum")
                        .long("freeze")
                        .help("Freeze the momentum field calculation")
                        .action(clap::ArgAction::SetTrue),
                ),
        )
        .subcommand(
            Command::new("post")
                .about("Post-process the simulation data")
                .subcommand(
                    Command::new("vtk")
                        .about("Write VTK files for visualization in Paraview")
                        .arg(
                            Arg::new("physical_data")
                                .long("physical")
                                .help("If the data is written in physical units")
                                .action(clap::ArgAction::SetTrue),
                        )
                        .arg(
                            Arg::new("node_type")
                                .long("node-type")
                                .help("Write node_type.vtk file")
                                .action(clap::ArgAction::SetTrue),
                        ),
                )
                .subcommand(
                    Command::new("unify")
                        .about(
                            "Unify the parallel output files into a single file for each time step",
                        )
                        .arg(
                            Arg::new("keep")
                                .short('k')
                                .long("keep")
                                .help("Keep the intermediate files after unification")
                                .action(clap::ArgAction::SetTrue),
                        ),
                ),
        )
        .get_matches();
    Ok(matches)
}

pub(crate) fn parse_matches(matches: &clap::ArgMatches) -> LbResult<Config> {
    let number_of_threads = *matches
        .get_one::<NonZeroUsize>("number_of_threads")
        .expect("Has 1 as default");
    let core_affinity = matches.get_flag("core_affinity");
    match matches.subcommand() {
        Some(("run", sub_m)) => {
            let frequency = *sub_m.get_one::<usize>("write_data").unwrap();
            let cfg = Config {
                mode: Mode::Run,
                number_of_threads,
                core_affinity,
                write_data: Some(WriteDataMode::Frequency(frequency)),
                max_iterations: sub_m.get_one::<usize>("max_iterations").cloned(),
                freeze_momentum: sub_m.get_flag("freeze_momentum"),
                ..Default::default()
            };
            Ok(cfg)
        }
        Some(("post", sub_m)) => match sub_m.subcommand() {
            Some(("vtk", sub_m)) => {
                let cfg = Config {
                    mode: Mode::PostVTK,
                    number_of_threads,
                    core_affinity,
                    node_type: sub_m.get_flag("node_type"),
                    physical_data: sub_m.get_flag("physical_data"),
                    ..Default::default()
                };
                Ok(cfg)
            }
            Some(("unify", sub_m)) => {
                let cfg = Config {
                    mode: Mode::PostUnify,
                    number_of_threads,
                    core_affinity,
                    keep: sub_m.get_flag("keep"),
                    ..Default::default()
                };
                Ok(cfg)
            }
            _ => unreachable!("At least one subcommand is required: .subcommand_required(true)"),
        },
        _ => unreachable!("At least one subcommand is required: .subcommand_required(true)"),
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
