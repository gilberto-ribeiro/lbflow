use crate::prelude::*;
use clap::{Arg, Command};
use core_affinity::{get_core_ids, set_for_current};
use std::num::{NonZero, NonZeroUsize};

pub type LbResult<T> = Result<T, Box<dyn std::error::Error>>;

#[derive(Debug)]
pub struct Config {
    pub mode: Mode,
    pub number_of_threads: NonZeroUsize,
    pub core_affinity: bool,
    pub write_data: Option<WriteDataMode>,
    pub max_iterations: Option<usize>,
    pub freeze_momentum: bool,
    pub node_type: bool,
    pub physical_data: bool,
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
        }
    }
}

impl Config {
    pub fn get_number_of_threads(&self) -> usize {
        usize::from(self.number_of_threads)
    }

    pub fn get_write_data_mode(&self) -> &WriteDataMode {
        self.write_data.as_ref().unwrap()
    }

    pub fn get_max_iterations(&self) -> usize {
        self.max_iterations.unwrap()
    }
}

#[derive(Debug)]
pub enum Mode {
    Run,
    Post,
}

pub fn get_args() -> LbResult<clap::ArgMatches> {
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
                .short('a')
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
                        .short('m')
                        .long("max-iterations")
                        .value_name("ITER")
                        .help("The maximum number of iterations")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("100000"),
                )
                .arg(
                    Arg::new("freeze_momentum")
                        .short('f')
                        .long("freeze")
                        .help("Freeze the momentum field calculation")
                        .action(clap::ArgAction::SetTrue),
                ),
        )
        .subcommand(
            Command::new("post")
                .about("Post-process the simulation data")
                .arg(
                    Arg::new("physical_data")
                        .short('p')
                        .long("physical")
                        .help("If the data is written in physical units")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("node_type")
                        .short('t')
                        .long("type")
                        .help("Write node_type.vtk file")
                        .action(clap::ArgAction::SetTrue),
                ),
        )
        .get_matches();
    Ok(matches)
}

pub fn parse_matches(matches: &clap::ArgMatches) -> LbResult<Config> {
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
        Some(("post", sub_m)) => {
            let cfg = Config {
                mode: Mode::Post,
                number_of_threads,
                core_affinity,
                node_type: sub_m.get_flag("node_type"),
                physical_data: sub_m.get_flag("physical_data"),
                ..Default::default()
            };
            Ok(cfg)
        }
        _ => unreachable!("At least one subcommand is required: .subcommand_required(true)"),
    }
}

pub fn init_global_pool(num_threads: usize, pin_all_cores: bool) {
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
