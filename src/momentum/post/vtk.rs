use crate::prelude::*;
use colored::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub fn read_coordinates_file(dim: usize) -> (Vec<Vec<usize>>, Vec<Vec<Float>>, Vec<NodeType>) {
    let path = Path::new(crate::io::DATA_PATH).join(crate::io::COORDINATES_FILE);
    let mut indexes: Vec<Vec<usize>> = Vec::new();
    let mut coordinates: Vec<Vec<Float>> = Vec::new();
    let mut node_types: Vec<NodeType> = Vec::new();
    crate::io::read_csv_file(path)
        .unwrap_or_else(|e| {
            eprintln!("Error: {e}");
            std::process::exit(1);
        })
        .iter()
        .for_each(|s| {
            let items = s.split(",").map(String::from).collect::<Vec<String>>();
            let index = items
                .get(0..dim)
                .unwrap()
                .iter()
                .map(|value| value.parse::<usize>().expect("Error"))
                .collect::<Vec<usize>>();
            let coordinate = items
                .get(dim..(2 * dim))
                .unwrap()
                .iter()
                .map(|value| value.parse::<Float>().expect("Error"))
                .collect::<Vec<Float>>();
            let node_type = match items.last().unwrap().as_str() {
                "fluid" => Fluid,
                "solid" => Solid,
                _ => panic!("Value not expected: must be 'solid' ou 'fluid'."),
            };
            indexes.push(index);
            coordinates.push(coordinate);
            node_types.push(node_type);
        });
    (indexes, coordinates, node_types)
}

fn read_densities(time_step: usize) -> Vec<Float> {
    let path = Path::new(crate::io::DATA_PATH)
        .join(time_step.to_string())
        .join(crate::io::DENSITY_FILE);
    crate::io::read_csv_file(path)
        .and_then(crate::io::parse_scalar_from_string)
        .unwrap_or_else(|e| {
            eprintln!("Error: {e}");
            std::process::exit(1);
        })
}

fn read_velocities(time_step: usize) -> Vec<Vec<Float>> {
    let path = Path::new(crate::io::DATA_PATH)
        .join(time_step.to_string())
        .join(crate::io::VELOCITY_FILE);
    crate::io::read_csv_file(path)
        .and_then(crate::io::parse_vector_from_string)
        .unwrap_or_else(|e| {
            eprintln!("Error: {e}");
            std::process::exit(1);
        })
}

fn unify_densities(time_step: usize, keep: bool) -> LbResult<()> {
    let dir = Path::new(crate::io::DATA_PATH).join(time_step.to_string());
    let prefix = "density";
    let header = "density";
    crate::io::unify_parallel_csv_files(dir, prefix, header, keep)
}

fn unify_velocities(time_step: usize, dim: usize, keep: bool) -> LbResult<()> {
    let dir = Path::new(crate::io::DATA_PATH).join(time_step.to_string());
    let prefix = "velocity";
    let directions = ["x", "y", "z"];
    let header = (0..dim)
        .map(|x| format!("velocity_{}", directions[x]))
        .collect::<Vec<String>>()
        .join(",");
    crate::io::unify_parallel_csv_files(dir, prefix, &header, keep)
}

fn write_node_type_vtk<P: AsRef<Path>>(
    path: P,
    n: &[usize],
    coordinates: &[Vec<Float>],
    node_types: &[NodeType],
) -> LbResult<()> {
    let point_data = n.iter().product::<usize>();
    let mut file = File::create(path)?;
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "lbflow simulation data")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET STRUCTURED_GRID")?;
    writeln!(
        file,
        "DIMENSIONS {} {} {}",
        n.first().unwrap_or(&1),
        n.get(1).unwrap_or(&1),
        n.get(2).unwrap_or(&1)
    )?;
    writeln!(file, "POINTS {point_data} float")?;
    coordinates.iter().for_each(|coordinate| {
        writeln!(
            file,
            "{:.4e} {:.4e} {:.4e}",
            coordinate[0],
            coordinate[1],
            coordinate.get(2).unwrap_or(&0.0)
        )
        .unwrap();
    });
    writeln!(file, "POINT_DATA {point_data}")?;
    writeln!(file, "SCALARS node_type float 1")?;
    writeln!(file, "LOOKUP_TABLE default")?;
    node_types
        .iter()
        .map(|node_type| match node_type {
            Solid => "1",
            Fluid => "0",
        })
        .for_each(|node_type| {
            writeln!(file, "{node_type}").unwrap();
        });
    Ok(())
}

fn write_momentum_vtk<P: AsRef<Path>>(
    path: P,
    config: &Config,
    conversion_factor: &momentum::ConversionFactor,
    n: &[usize],
    coordinates: &[Vec<Float>],
    densities: &[Float],
    velocities: &[Vec<Float>],
) -> LbResult<()> {
    let point_data = n.iter().product::<usize>();
    let mut file = File::create(path)?;
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "lbflow simulation data")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET STRUCTURED_GRID")?;
    writeln!(
        file,
        "DIMENSIONS {} {} {}",
        n.first().unwrap_or(&1),
        n.get(1).unwrap_or(&1),
        n.get(2).unwrap_or(&1)
    )?;
    writeln!(file, "POINTS {point_data} float")?;
    coordinates.iter().for_each(|coordinate| {
        writeln!(
            file,
            "{:.4e} {:.4e} {:.4e}",
            coordinate[0],
            coordinate[1],
            coordinate.get(2).unwrap_or(&0.0)
        )
        .unwrap();
    });
    writeln!(file, "POINT_DATA {point_data}")?;
    writeln!(file, "SCALARS density float 1")?;
    writeln!(file, "LOOKUP_TABLE default")?;
    densities.iter().for_each(|density| {
        writeln!(file, "{density:.8e}").unwrap();
    });
    writeln!(file, "VECTORS velocity float")?;
    velocities.iter().for_each(|velocity| {
        writeln!(
            file,
            "{:.8e} {:.8e} {:.8e}",
            velocity[0],
            velocity[1],
            velocity.get(2).unwrap_or(&0.0)
        )
        .unwrap();
    });
    if config.physical_data {
        writeln!(file, "SCALARS physical_pressure float 1")?;
        writeln!(file, "LOOKUP_TABLE default")?;
        densities.iter().for_each(|&density| {
            let physical_pressure = compute_physical_pressure(density, conversion_factor);
            writeln!(file, "{physical_pressure:.8e}").unwrap();
        });
        writeln!(file, "VECTORS physical_velocity float")?;
        velocities.iter().for_each(|velocity| {
            let physical_velocity = compute_physical_velocity(velocity, conversion_factor);
            writeln!(
                file,
                "{:.8e} {:.8e} {:.8e}",
                physical_velocity[0],
                physical_velocity[1],
                physical_velocity.get(2).unwrap_or(&0.0)
            )
            .unwrap();
        });
    }
    Ok(())
}

pub fn post_vtk(config: Config, momentum_params: momentum::Parameters) {
    let n = momentum_params.n.clone();
    let dim = n.len();
    let (_, coordinates, node_types) = read_coordinates_file(dim);
    let conversion_factor = momentum::ConversionFactor::from(momentum_params);
    node_type_vtk(&config, &n, &coordinates, &node_types);
    momentum_vtk(&config, &conversion_factor, &n, &coordinates);
}

pub fn post_unify(config: Config, momentum_params: momentum::Parameters) {
    let keep = config.keep;
    let dim = momentum_params.n.len();
    momentum_unify(dim, keep);
}

fn compute_physical_pressure(
    density: Float,
    conversion_factor: &momentum::ConversionFactor,
) -> Float {
    let density_prime = density - LATTICE_DENSITY;
    let pressure_prime = CS_2 * density_prime;
    conversion_factor.reference_pressure
        + pressure_prime * conversion_factor.pressure_conversion_factor
}

fn compute_physical_velocity(
    velocity: &[Float],
    conversion_factor: &momentum::ConversionFactor,
) -> Vec<Float> {
    velocity
        .iter()
        .map(|u_x| u_x * conversion_factor.velocity_conversion_factor)
        .collect::<Vec<Float>>()
}

pub fn momentum_vtk(
    config: &Config,
    conversion_factor: &momentum::ConversionFactor,
    n: &[usize],
    coordinates: &[Vec<Float>],
) {
    crate::io::collect_time_steps()
        .unwrap_or_else(|e| {
            eprintln!("Error: {e}");
            std::process::exit(1);
        })
        .into_par_iter()
        .for_each(|time_step| {
            let case_name = crate::io::get_case_name();
            let file_name = format!("{case_name}_output_{time_step:08}.vtk");
            let path = Path::new(crate::io::VTK_PATH).join(&file_name);
            let densities = read_densities(time_step);
            let velocities = read_velocities(time_step);
            println!(
                "Writing {} for time step {}.",
                file_name.bold().yellow(),
                time_step.to_string().bold().yellow()
            );
            if let Err(e) = write_momentum_vtk(
                &path,
                config,
                conversion_factor,
                n,
                coordinates,
                &densities,
                &velocities,
            ) {
                eprintln!("Error writing VTK {path:?}: {e}");
                std::process::exit(1);
            };
        });
}

pub fn momentum_unify(dim: usize, keep: bool) {
    crate::io::collect_time_steps()
        .unwrap_or_else(|e| {
            eprintln!("Error: {e}");
            std::process::exit(1);
        })
        .into_par_iter()
        .for_each(|time_step| {
            println!(
                "Unifying {} and {} files for time step {}.",
                "density".bold().yellow(),
                "velocity".bold().yellow(),
                time_step.to_string().bold().yellow()
            );
            if let Err(e) = unify_densities(time_step, keep) {
                eprintln!("Error unifying densities for time step {time_step}: {e}");
                std::process::exit(1);
            };
            if let Err(e) = unify_velocities(time_step, dim, keep) {
                eprintln!("Error unifying velocities for time step {time_step}: {e}");
                std::process::exit(1);
            };
        });
}

pub fn node_type_vtk(
    config: &Config,
    n: &[usize],
    coordinates: &[Vec<Float>],
    node_types: &[NodeType],
) {
    if config.node_type {
        let path = Path::new(crate::io::VTK_PATH).join(crate::io::NODE_TYPE_VTK_FILE);
        println!(
            "Writing {} for node types.\n",
            crate::io::NODE_TYPE_VTK_FILE.bold().yellow()
        );
        if let Err(e) = write_node_type_vtk(&path, n, coordinates, node_types) {
            eprintln!("Error writing VTK {path:?}: {e}");
            std::process::exit(1);
        }
    }
}
