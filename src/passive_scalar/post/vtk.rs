use crate::prelude::*;
use colored::*;
use rayon::prelude::*;
use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;

fn read_concentrations(time_step: usize, prefix: &str) -> Vec<Float> {
    let dir = Path::new(crate::io::DATA_PATH).join(time_step.to_string());
    crate::io::read_parallel_csv_files(dir, prefix)
        .and_then(crate::io::parse_scalar_from_string)
        .unwrap_or_else(|e| {
            eprintln!("Error: {e}");
            std::process::exit(1);
        })
}

fn append_passive_scalar_vtk<P: AsRef<Path>>(
    path: P,
    prefix: &str,
    _config: &Config,
    concentrations: &[Float],
) -> LbResult<()> {
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    writeln!(file, "SCALARS {prefix} float 1")?;
    writeln!(file, "LOOKUP_TABLE default")?;
    concentrations.iter().for_each(|concentration| {
        writeln!(file, "{concentration:.8e}").unwrap();
    });
    Ok(())
}

fn passive_scalar_vtk(config: &Config, prefix: &str) {
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
            let concentrations = read_concentrations(time_step, prefix);
            println!(
                "Appending {} in {} for time step {}.",
                prefix.bold().yellow(),
                file_name.bold().yellow(),
                time_step.to_string().bold().yellow()
            );
            if let Err(e) = append_passive_scalar_vtk(&path, prefix, config, &concentrations) {
                eprintln!("Error writing VTK {path:?}: {e}");
                std::process::exit(1);
            };
        });
}

pub fn post_vtk(
    config: Config,
    momentum_params: momentum::Parameters,
    passive_scalar_params: passive_scalar::Parameters,
) {
    let n = momentum_params.n.clone();
    let dim = n.len();
    let (_, coordinates, node_types) = momentum::post::vtk::read_coordinates_file(dim);
    let conversion_factor = momentum::ConversionFactor::from(momentum_params);
    momentum::post::vtk::node_type_vtk(&config, &n, &coordinates, &node_types);
    momentum::post::vtk::momentum_vtk(&config, &conversion_factor, &n, &coordinates);
    passive_scalar_vtk(&config, &passive_scalar_params.scalar_name);
}
