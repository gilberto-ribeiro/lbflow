use crate::prelude_crate::*;
use colored::*;
use rayon::prelude::*;
use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;

fn read_scalar_values(time_step: usize, prefix: &str) -> Vec<Float> {
    let file_name = format!("{prefix}.csv");
    let path = Path::new(crate::io::DATA_PATH)
        .join(time_step.to_string())
        .join(file_name);
    crate::io::read_csv_file(path)
        .and_then(crate::io::parse_scalar_from_string)
        .unwrap_or_else(|e| {
            eprintln!("Error: {e}");
            std::process::exit(1);
        })
}

fn unify_scalar_values(time_step: usize, prefix: &str, keep: bool) -> LbResult<()> {
    let dir = Path::new(crate::io::DATA_PATH).join(time_step.to_string());
    let header = prefix;
    crate::io::unify_parallel_csv_files(dir, prefix, header, keep)
}

fn append_passive_scalar_vtk<P: AsRef<Path>>(
    path: P,
    prefix: &str,
    scalar_values: &[Float],
) -> LbResult<()> {
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    writeln!(file, "SCALARS {prefix} float 1")?;
    writeln!(file, "LOOKUP_TABLE default")?;
    scalar_values.iter().for_each(|scalar_value| {
        writeln!(file, "{scalar_value:.8e}").unwrap();
    });
    Ok(())
}

fn passive_scalar_vtk(prefix: &str) {
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
            let scalar_values = read_scalar_values(time_step, prefix);
            println!(
                "Appending {} in {} for time step {}.",
                prefix.bold().yellow(),
                file_name.bold().yellow(),
                time_step.to_string().bold().yellow()
            );
            if let Err(e) = append_passive_scalar_vtk(&path, prefix, &scalar_values) {
                eprintln!("Error writing VTK {path:?}: {e}");
                std::process::exit(1);
            };
        });
}

fn passive_scalar_unify(prefix: &str, keep: bool) {
    crate::io::collect_time_steps()
        .unwrap_or_else(|e| {
            eprintln!("Error: {e}");
            std::process::exit(1);
        })
        .into_par_iter()
        .for_each(|time_step| {
            println!(
                "Unifying {} files for time step {}.",
                prefix.bold().yellow(),
                time_step.to_string().bold().yellow()
            );
            if let Err(e) = unify_scalar_values(time_step, prefix, keep) {
                eprintln!("Error unifying {prefix} values for time step {time_step}: {e}");
                std::process::exit(1);
            };
        });
}

fn passive_scalar_vtk_vec(prefixes: &[String]) {
    prefixes.iter().for_each(|prefix| {
        passive_scalar_vtk(prefix);
    });
}

fn passive_scalar_unify_vec(prefixes: &[String], keep: bool) {
    prefixes.iter().for_each(|prefix| {
        passive_scalar_unify(prefix, keep);
    });
}

pub(crate) fn post_vtk(
    cli_args: Cli,
    momentum_params: momentum::Parameters,
    passive_scalar_params: passive_scalar::Parameters,
) {
    match cli_args.command {
        cli::Command::Run { .. } => {}
        cli::Command::Post { command } => {
            if let cli::PostCommand::Vtk {
                node_type,
                physical_data,
                ..
            } = command
            {
                let n = momentum_params.n.clone();
                let dim = n.len();
                let (_, coordinates, node_types) = momentum::post::vtk::read_coordinates_file(dim);
                let conversion_factor = momentum::ConversionFactor::from(momentum_params);
                momentum::post::vtk::node_type_vtk(&n, &coordinates, &node_types, node_type);
                momentum::post::vtk::momentum_vtk(
                    &conversion_factor,
                    &n,
                    &coordinates,
                    physical_data,
                );
                passive_scalar_vtk(passive_scalar_params.scalar_name);
            }
        }
    }
}

pub(crate) fn post_unify(
    cli_args: Cli,
    momentum_params: momentum::Parameters,
    passive_scalar_params: passive_scalar::Parameters,
) {
    match cli_args.command {
        cli::Command::Run { .. } => {}
        cli::Command::Post { command } => {
            if let cli::PostCommand::Unify { keep } = command {
                let dim = momentum_params.n.len();
                momentum::post::vtk::momentum_unify(dim, keep);
                passive_scalar_unify(passive_scalar_params.scalar_name, keep);
            }
        }
    }
}

pub(crate) fn post_vtk_vec(
    cli_args: Cli,
    momentum_params: momentum::Parameters,
    passive_scalar_params_vec: Vec<passive_scalar::Parameters>,
) {
    match cli_args.command {
        cli::Command::Run { .. } => {}
        cli::Command::Post { command } => {
            if let cli::PostCommand::Vtk {
                node_type,
                physical_data,
                ..
            } = command
            {
                let n = momentum_params.n.clone();
                let dim = n.len();
                let (_, coordinates, node_types) = momentum::post::vtk::read_coordinates_file(dim);
                let conversion_factor = momentum::ConversionFactor::from(momentum_params);
                momentum::post::vtk::node_type_vtk(&n, &coordinates, &node_types, node_type);
                momentum::post::vtk::momentum_vtk(
                    &conversion_factor,
                    &n,
                    &coordinates,
                    physical_data,
                );
                passive_scalar_vtk_vec(
                    &passive_scalar_params_vec
                        .iter()
                        .map(|passive_scalar_params| passive_scalar_params.scalar_name.to_string())
                        .collect::<Vec<String>>(),
                );
            }
        }
    }
}

pub(crate) fn post_unify_vec(
    cli_args: Cli,
    momentum_params: momentum::Parameters,
    passive_scalar_params_vec: Vec<passive_scalar::Parameters>,
) {
    match cli_args.command {
        cli::Command::Run { .. } => {}
        cli::Command::Post { command } => {
            if let cli::PostCommand::Unify { keep } = command {
                let dim = momentum_params.n.len();
                momentum::post::vtk::momentum_unify(dim, keep);
                passive_scalar_unify_vec(
                    &passive_scalar_params_vec
                        .iter()
                        .map(|passive_scalar_params| passive_scalar_params.scalar_name.to_string())
                        .collect::<Vec<String>>(),
                    keep,
                );
            }
        }
    }
}
