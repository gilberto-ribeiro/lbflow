pub mod vtk;

use crate::prelude::*;
use rayon::prelude::*;

pub type PostComputation = fn(&momentum::Lattice) -> Vec<PostResult>;

pub struct PostResult {
    pub name: String,
    pub label: String,
    pub value: Float,
    pub unit: Option<String>,
}

impl PostResult {
    pub fn new(name: String, label: String, value: Float, unit: Option<String>) -> Self {
        Self {
            name,
            label,
            value,
            unit,
        }
    }
}

#[derive(Debug)]
pub struct PostFunction {
    pub file_name: String,
    pub interval: usize,
    pub function: PostComputation,
}

impl PostFunction {
    pub fn new(file_name: String, interval: usize, function: PostComputation) -> Self {
        Self {
            file_name,
            interval,
            function,
        }
    }
}

pub fn compute_mean_density(lattice: &momentum::Lattice) -> Vec<PostResult> {
    let rho_sum = lattice
        .get_fluid_nodes()
        .par_iter()
        .map(|node| node.get_density())
        .sum::<Float>();
    let number_of_fluid_nodes = lattice.get_fluid_nodes().len() as Float;
    let rho_mean = rho_sum / number_of_fluid_nodes;
    let rho_result: PostResult = PostResult::new(
        "mean_density".to_string(),
        "mean density".to_string(),
        rho_mean,
        None,
    );
    vec![rho_result]
}
