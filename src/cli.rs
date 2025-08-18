use crate::prelude::*;
use clap::{Arg, Command};
use std::num::{NonZero, NonZeroUsize};

pub type LbResult<T> = Result<T, Box<dyn std::error::Error>>;

#[derive(Debug)]
pub struct Config {
    pub mode: Mode,
    pub number_of_threads: NonZeroUsize,
    pub case_name: Option<String>,
    pub write_data: Option<WriteDataMode>,
    pub max_iterations: Option<usize>,
    pub node_type: bool,
    pub physical_data: bool,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            mode: Mode::Run,
            number_of_threads: NonZero::new(1).unwrap(),
            case_name: None,
            write_data: None,
            max_iterations: None,
            node_type: false,
            physical_data: false,
        }
    }
}

impl Config {
    pub fn get_number_of_threads(&self) -> usize {
        usize::from(self.number_of_threads)
    }

    pub fn get_case_name(&self) -> &str {
        self.case_name.as_deref().unwrap_or("unknown")
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
        .subcommand(
            Command::new("run")
                .about("Run the simulation")
                .arg(
                    Arg::new("case_name")
                        .short('c')
                        .long("case-name")
                        .value_name("CASE")
                        .help("The name of the case to be simulated")
                        .required(true),
                )
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
    match matches.subcommand() {
        Some(("run", sub_m)) => {
            let frequency = *sub_m.get_one::<usize>("write_data").unwrap();
            let cfg = Config {
                mode: Mode::Run,
                number_of_threads,
                case_name: sub_m.get_one::<String>("case_name").cloned(),
                write_data: Some(WriteDataMode::Frequency(frequency)),
                max_iterations: sub_m.get_one::<usize>("max_iterations").cloned(),
                ..Default::default()
            };
            Ok(cfg)
        }
        Some(("post", sub_m)) => {
            let cfg = Config {
                mode: Mode::Post,
                number_of_threads,
                node_type: sub_m.get_flag("node_type"),
                physical_data: sub_m.get_flag("physical_data"),
                ..Default::default()
            };
            Ok(cfg)
        }
        _ => unreachable!("At least one subcommand is required: .subcommand_required(true)"),
    }
}
