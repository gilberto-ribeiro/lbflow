// ------------------------------------------------------------------------------- IMPORTS

use crate::prelude_crate::*;
use colored::*;
use core::panic;
use regex::Regex;
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process;

pub(crate) const DATA_PATH: &str = "./data";
pub(crate) const PRE_PROCESSING_PATH: &str = "./pre_processing";
pub(crate) const POST_PROCESSING_PATH: &str = "./post_processing";
pub(crate) const VTK_PATH: &str = "./post_processing/vtk_files";
pub(crate) const COORDINATES_FILE: &str = "coordinates.csv";
pub(crate) const DENSITY_FILE: &str = "density.csv";
pub(crate) const VELOCITY_FILE: &str = "velocity.csv";
pub(crate) const RESIDUALS_FILE: &str = "residuals.csv";
pub(crate) const NODE_TYPE_VTK_FILE: &str = "node_type.vtk";
pub(crate) const RESIDUALS_GRAPH_FILE: &str = "gr_residuals.gp";
pub(crate) const LIVE_RESIDUALS_GRAPH_FILE: &str = "live_residuals.gp";
pub(crate) const BOUNCE_BACK_MAP_FILE: &str = "map.xyz";

pub struct ResidualsInfo {
    pub(crate) print_header: String,
    pub(crate) print_line: String,
    pub(crate) write_header: String,
    pub(crate) write_line: String,
    pub(crate) time_step: usize,
}

#[derive(Debug)]
pub enum _WriteDataMode {
    Frequency(usize),
    ListOfSteps(Vec<usize>),
}

pub(crate) fn create_case_directories() -> LbResult<()> {
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

pub(crate) fn read_bounce_back_map() -> Vec<NodeType> {
    let pre_processing_path = Path::new(crate::io::PRE_PROCESSING_PATH);
    let path = pre_processing_path.join(crate::io::BOUNCE_BACK_MAP_FILE);
    let data = read_xyz_file(path).unwrap_or_else(|_| {
        eprintln!("Error while reading the bounce back map file.");
        process::exit(1);
    });
    parse_node_type_from_string(&data)
}

pub(crate) fn progress_bar(current: usize, total: usize) {
    let current = current + 1;
    let percentage = current as Float / total as Float;
    let bar_length = 50;
    let filled_length = (bar_length as Float * percentage) as usize;
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

fn parse_node_type_from_string(data: &[String]) -> Vec<NodeType> {
    data.iter()
        .map(|s| match s.as_str() {
            "0" => NodeType::Fluid,
            "1" => NodeType::Solid,
            _ => panic!("Invalid node type: {s}"),
        })
        .collect()
}

pub(crate) fn parse_scalar_from_string(data: Vec<String>) -> LbResult<Vec<Float>> {
    data.iter()
        .map(|s| {
            s.parse::<Float>()
                .map_err(|e| format!("Error parsing value `{s}`: {e}").into())
        })
        .collect()
}

pub(crate) fn parse_vector_from_string(data: Vec<String>) -> LbResult<Vec<Vec<Float>>> {
    data.iter()
        .map(|s| {
            s.split(',')
                .map(|value| {
                    value
                        .parse::<Float>()
                        .map_err(|e| format!("Error parsing value `{value}`: {e}").into())
                })
                .collect()
        })
        .collect()
}

fn read_xyz_file<P: AsRef<Path>>(path: P) -> LbResult<Vec<String>> {
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    let vector = contents
        .trim()
        .split_terminator("\n\n")
        .flat_map(|slice| {
            slice
                .trim()
                .lines()
                .map(|line| {
                    line.split_whitespace()
                        .map(|value| value.to_string())
                        .collect::<Vec<String>>()
                })
                .rev()
                .collect::<Vec<Vec<String>>>()
                .concat()
        })
        .collect::<Vec<String>>();
    Ok(vector)
}

pub(crate) fn read_csv_file<P: AsRef<Path>>(path: P) -> LbResult<Vec<String>> {
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    let mut contents = contents.trim().split_terminator('\n');
    contents.next();
    let data = contents
        .map(|line| line.to_string())
        .collect::<Vec<String>>();
    Ok(data)
}

pub(crate) fn collect_time_steps() -> LbResult<Vec<usize>> {
    let mut time_steps = Vec::new();
    let dir = DATA_PATH;
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        if entry.file_type()?.is_dir() {
            if let Some(name) = entry.file_name().to_str() {
                if let Ok(time_step) = name.parse::<usize>() {
                    time_steps.push(time_step);
                }
            }
        }
    }
    time_steps.sort_unstable();
    Ok(time_steps)
}

fn collect_parallel_files<P: AsRef<Path>>(dir: P, prefix: &str) -> LbResult<Vec<PathBuf>> {
    let pattern = format!(r"^{}_(\d+)\.csv$", regex::escape(prefix));
    let re = Regex::new(&pattern).unwrap();
    let mut items: Vec<(usize, PathBuf)> = Vec::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let file_name = entry.file_name();
        let file_name = match file_name.to_str() {
            Some(s) => s,
            None => continue,
        };
        if let Some(caps) = re.captures(file_name) {
            let idx: usize = caps[1].parse().unwrap();
            items.push((idx, entry.path()));
        }
    }
    items.sort_by_key(|(idx, _)| *idx);
    Ok(items.into_iter().map(|(_, p)| p).collect())
}

fn delete_parallel_files<P: AsRef<Path>>(dir: P, prefix: &str) -> LbResult<()> {
    let paths = collect_parallel_files(dir, prefix)?;
    for path in paths {
        fs::remove_file(path)?;
    }
    Ok(())
}

fn read_parallel_csv_files<P: AsRef<Path>>(dir: P, prefix: &str) -> LbResult<Vec<String>> {
    Ok(collect_parallel_files(dir, prefix)?
        .into_iter()
        .flat_map(|path| match read_csv_file(path) {
            Ok(data) => data,
            Err(e) => {
                eprintln!("Error reading CSV file: {e}");
                process::exit(1);
            }
        })
        .collect())
}

pub(crate) fn unify_parallel_csv_files<P: AsRef<Path>>(
    dir: P,
    prefix: &str,
    header: &str,
    keep: bool,
) -> LbResult<()> {
    let file_name = format!("{prefix}.csv");
    let path = dir.as_ref().join(&file_name);
    let data = read_parallel_csv_files(dir.as_ref(), prefix)?;
    let mut file = File::create(&path)?;
    writeln!(file, "{header}")?;
    data.iter()
        .for_each(|line| writeln!(file, "{line}").unwrap());
    if !keep {
        delete_parallel_files(dir.as_ref(), prefix)?;
    }
    Ok(())
}

pub(crate) fn print_residuals(info: &ResidualsInfo) {
    if info.time_step.is_multiple_of(100) {
        println!("{}\n", info.print_header);
    }
    println!("{}", info.print_line);
}

pub(crate) fn write_residuals(info: &ResidualsInfo) -> LbResult<()> {
    let data_path = Path::new(crate::io::DATA_PATH);
    let path = data_path.join(crate::io::RESIDUALS_FILE);
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    if info.time_step == 0 {
        writeln!(file, "{}", info.write_header)?;
    }
    writeln!(file, "{}", info.write_line)?;
    Ok(())
}

pub(crate) fn get_case_name() -> String {
    let exe = std::env::current_exe().unwrap();
    Path::new(&exe)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(env!("CARGO_MANIFEST_DIR"))
        .to_string()
}

pub(crate) fn generate_uniform_scalars(value: Float, n: &[usize]) -> Vec<Float> {
    let num_nodes = n.iter().product();
    vec![value; num_nodes]
}

pub(crate) fn generate_uniform_vectors(value: Vec<Float>, n: &[usize]) -> Vec<Vec<Float>> {
    let num_nodes = n.iter().product();
    vec![value; num_nodes]
}

pub(crate) fn generate_scalars_from_file<P>(path: P) -> Vec<Float>
where
    P: AsRef<Path>,
{
    let error_msg = format!("Error reading {}.", path.as_ref().display());
    let data = crate::io::read_csv_file(path).unwrap_or_else(|_| {
        eprintln!("{}", error_msg);
        std::process::exit(1);
    });
    crate::io::parse_scalar_from_string(data).unwrap()
}

pub(crate) fn generate_vectors_from_file<P>(path: P) -> Vec<Vec<Float>>
where
    P: AsRef<Path>,
{
    let error_msg = format!("Error reading {}.", path.as_ref().display());
    let data = crate::io::read_csv_file(path).unwrap_or_else(|_| {
        eprintln!("{}", error_msg);
        std::process::exit(1);
    });
    crate::io::parse_vector_from_string(data).unwrap()
}

pub(crate) fn generate_scalars_from_time_step(time_step: usize, filename: &str) -> Vec<Float> {
    let data_path = Path::new(crate::io::DATA_PATH);
    let step_path = data_path.join(time_step.to_string());
    let path = step_path.join(filename);
    generate_scalars_from_file(path)
}

pub(crate) fn generate_vectors_from_time_step(time_step: usize, filename: &str) -> Vec<Vec<Float>> {
    let data_path = Path::new(crate::io::DATA_PATH);
    let step_path = data_path.join(time_step.to_string());
    let path = step_path.join(filename);
    generate_vectors_from_file(path)
}
