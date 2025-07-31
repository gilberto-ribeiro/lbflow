use crate::NodeType;
use crate::momentum::Node;
use colored::*;
use core::panic;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::Path;
use std::process;

pub const DATA_PATH: &str = "./data";
pub const PRE_PROCESSING_PATH: &str = "./pre_processing";
pub const POST_PROCESSING_PATH: &str = "./post_processing";
pub const VTK_PATH: &str = "./post_processing/vtk_files";
pub const COORDINATES_FILE: &str = "coordinates.csv";
pub const DENSITY_FILE: &str = "density.csv";
pub const VELOCITY_FILE: &str = "velocity.csv";
pub const RESIDUALS_FILE: &str = "residuals.csv";
pub const RESIDUALS_GRAPH_FILE: &str = "gr_residuals.gp";
pub const LIVE_RESIDUALS_GRAPH_FILE: &str = "live_residuals.gp";
pub const BOUNCE_BACK_MAP_FILE: &'static str = "map.dat";

pub enum WriteDataMode {
    Frequency(usize),
    ListOfSteps(Vec<usize>),
}

pub fn create_case_directories() -> io::Result<()> {
    let list_of_paths = [
        DATA_PATH,
        PRE_PROCESSING_PATH,
        POST_PROCESSING_PATH,
        VTK_PATH,
    ];
    for path_str in list_of_paths {
        let path = Path::new(path_str);
        if !path.exists() {
            println!("Creating the {} path.\n", path_str.yellow().bold());
            if let Err(e) = fs::create_dir(path) {
                eprintln!("Error while creating the {path_str} path: {e}.");
                process::exit(1);
            };
        } else {
            println!("The {} path already exists.\n", path_str.yellow().bold());
        }
    }
    Ok(())
}

pub fn read_bounce_back_map() -> io::Result<Vec<NodeType>> {
    let pre_processing_path = Path::new(crate::io::PRE_PROCESSING_PATH);
    let path = pre_processing_path.join(crate::io::BOUNCE_BACK_MAP_FILE);
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    let map = contents
        .trim()
        .split_terminator("\n\n")
        .map(|slice| {
            let mut slice = slice
                .trim()
                .lines()
                .map(|line| {
                    line.split_whitespace()
                        .map(|value| {
                            let value = value.parse::<i32>().unwrap();
                            match value {
                                0 => NodeType::Fluid,
                                1 => NodeType::Solid,
                                _ => panic!("Invalid value in bounce-back map: {value}"),
                            }
                        })
                        .collect::<Vec<NodeType>>()
                })
                .collect::<Vec<Vec<NodeType>>>();
            slice.reverse();
            slice.concat()
        })
        .collect::<Vec<Vec<NodeType>>>()
        .concat();
    Ok(map)
}

pub fn progress_bar(current: usize, total: usize) {
    let current = current + 1;
    let percentage = current as f64 / total as f64;
    let bar_length = 50;
    let filled_length = (bar_length as f64 * percentage) as usize;
    let completed = "█".repeat(filled_length);
    let remaining = "░".repeat(bar_length - filled_length);
    let bar = completed + &remaining;
    let percentage = percentage * 100.0;
    print!("\r{}", format!("{bar} {percentage:.2}%").green().bold());
    std::io::Write::flush(&mut std::io::stdout()).unwrap();
    if current == total {
        println!();
        println!();
    }
}
