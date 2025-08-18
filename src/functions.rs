use crate::prelude::*;
use std::path::Path;

pub fn from_bounce_back_map_file() -> Vec<NodeType> {
    crate::io::read_bounce_back_map()
}

pub fn only_fluid_nodes(n: Vec<usize>) -> Vec<NodeType> {
    let num_nodes = n.iter().product();
    vec![NodeType::Fluid; num_nodes]
}

pub fn density_from_time_step(time_step: usize) -> Vec<Float> {
    let data_path = Path::new(crate::io::DATA_PATH);
    let step_path = data_path.join(time_step.to_string());
    let path = step_path.join(crate::io::DENSITY_FILE);
    from_density_file(path)
}

pub fn velocity_from_time_step(time_step: usize) -> Vec<Vec<Float>> {
    let data_path = Path::new(crate::io::DATA_PATH);
    let step_path = data_path.join(time_step.to_string());
    let path = step_path.join(crate::io::VELOCITY_FILE);
    from_velocity_file(path)
}

pub fn from_density_file<P>(path: P) -> Vec<Float>
where
    P: AsRef<Path>,
{
    let data = crate::io::read_csv_file(path).unwrap_or_else(|_| {
        eprintln!("Error reading the density file.");
        std::process::exit(1);
    });
    crate::io::parse_scalar_from_string(data).unwrap()
}

pub fn from_velocity_file<P>(path: P) -> Vec<Vec<Float>>
where
    P: AsRef<Path>,
{
    let data = crate::io::read_csv_file(path).unwrap_or_else(|_| {
        eprintln!("Error reading the velocity file.");
        std::process::exit(1);
    });
    crate::io::parse_vector_from_string(data).unwrap()
}

pub fn uniform_density(value: Float, n: Vec<usize>) -> Vec<Float> {
    uniform_scalar(value, n)
}

pub fn uniform_velocity(value: Vec<Float>, n: Vec<usize>) -> Vec<Vec<Float>> {
    uniform_vector(value, n)
}

fn uniform_scalar(value: Float, n: Vec<usize>) -> Vec<Float> {
    let num_nodes = n.iter().product();
    vec![value; num_nodes]
}

fn uniform_vector(value: Vec<Float>, n: Vec<usize>) -> Vec<Vec<Float>> {
    let num_nodes = n.iter().product();
    vec![value; num_nodes]
}
