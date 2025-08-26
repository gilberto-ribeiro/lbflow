// ------------------------------------------------------------------------------- IMPORTS

use super::Lattice;
use crate::prelude::*;
use colored::*;
use rayon::prelude::*;
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::process;

impl Lattice {
    pub fn get_residuals_info(&self) -> crate::io::ResidualsInfo {
        crate::io::ResidualsInfo {
            print_header: self.print_residuals_header(),
            print_line: self.print_residuals_line(),
            write_header: self.write_residuals_header(),
            write_line: self.write_residuals_line(),
            time_step: self.get_time_step(),
        }
    }

    fn print_residuals_header(&self) -> String {
        let binding = self.get_residuals();
        let residuals_velocity = binding.get_velocity();
        let directions = ["x", "y", "z"];
        let mut print_header = format!(
            "\n{:>10} {:>16}",
            "time_step".cyan().bold(),
            "density".cyan().bold()
        );
        residuals_velocity
            .iter()
            .zip(directions.iter())
            .for_each(|(_, x)| {
                print_header.push_str(&format!(" {:>16}", format!("velocity_{x}").cyan().bold()));
            });
        print_header
    }

    fn print_residuals_line(&self) -> String {
        let residuals_density = self.get_residuals().get_density();
        let binding = self.get_residuals();
        let residuals_velocity = binding.get_velocity();
        let mut print_line = format! {"{:>10} {:>16.8e}", self.get_time_step(), residuals_density};
        residuals_velocity.iter().for_each(|u_x| {
            print_line.push_str(&format!(" {u_x:>16.8e}"));
        });
        print_line
    }

    fn write_residuals_header(&self) -> String {
        let mut write_header = String::from("time_step,density");
        let directions = ["x", "y", "z"];
        let velocities = (0..*self.get_d())
            .map(|x| format!("velocity_{}", directions[x]))
            .collect::<Vec<String>>()
            .join(",");
        write_header = write_header + "," + &velocities;
        write_header
    }

    fn write_residuals_line(&self) -> String {
        let residuals_density = self.get_residuals().get_density();
        let binding = self.get_residuals();
        let residuals_velocity = binding.get_velocity();
        let mut write_line = format! {"{},{:.8e}", self.get_time_step(), residuals_density};
        residuals_velocity.iter().for_each(|u_x| {
            write_line.push_str(&format!(",{u_x:.8e}"));
        });
        write_line
    }

    pub fn write_data(&self) {
        match self.get_config().get_write_data_mode() {
            WriteDataMode::Frequency(n) => {
                if self.get_time_step() % n == 0 || self.get_time_step() == 0 {
                    println!();
                    self.write_data_from_steps();
                }
            }
            WriteDataMode::ListOfSteps(list) => {
                if list.contains(&self.get_time_step()) || self.get_time_step() == 0 {
                    println!();
                    self.write_data_from_steps();
                }
            }
        }
    }

    fn write_data_from_steps(&self) {
        let data_path = Path::new(crate::io::DATA_PATH);
        let step_path = data_path.join(self.get_time_step().to_string());
        if let Err(e) = fs::create_dir_all(&step_path) {
            eprintln!(
                "Error while creating the the step {} directory: {e}.",
                self.get_time_step()
            );
            process::exit(1);
        };
        println!(
            "\nWriting {} for time step {}.\n",
            crate::io::DENSITY_FILE.yellow().bold(),
            self.get_time_step().to_string().yellow().bold()
        );
        if let Err(e) = self.write_density(&step_path) {
            eprintln!("Error while writing the density file: {e}.");
            process::exit(1);
        };
        println!(
            "\nWriting {} for time step {}.\n",
            crate::io::VELOCITY_FILE.yellow().bold(),
            self.get_time_step().to_string().yellow().bold(),
        );
        if let Err(e) = self.write_velocity(&step_path) {
            eprintln!("Error while writing the velocity file: {e}.");
            process::exit(1);
        };
    }

    fn write_density<P>(&self, step_path: P) -> LbResult<()>
    where
        P: AsRef<Path>,
    {
        let number_of_threads = self.get_config().get_number_of_threads();
        let step_path = Arc::new(step_path.as_ref().to_path_buf());
        if number_of_threads == 1 {
            let path = step_path.join(crate::io::DENSITY_FILE);
            let mut file = File::create(path)?;
            writeln!(file, "density")?;
            self.get_nodes().iter().for_each(|node| {
                let line = format!("{:.8e}", node.get_density());
                writeln!(file, "{line}").unwrap();
            });
        } else if number_of_threads > 1 {
            let number_of_nodes = self.get_number_of_nodes();
            let chunk_size = number_of_nodes / number_of_threads;
            self.get_nodes()
                .par_chunks(chunk_size)
                .enumerate()
                .for_each(|(i, chunk)| {
                    let file_str = format!("density_{i}.csv");
                    let path = step_path.join(file_str);
                    let mut file = File::create(path).unwrap();
                    writeln!(file, "density").unwrap();
                    chunk.iter().for_each(|node| {
                        let line = format!("{:.8e}", node.get_density());
                        writeln!(file, "{line}").unwrap();
                    });
                });
        }
        Ok(())
    }

    fn write_velocity<P>(&self, step_path: P) -> LbResult<()>
    where
        P: AsRef<Path>,
    {
        let number_of_threads = self.get_config().get_number_of_threads();
        let step_path = Arc::new(step_path.as_ref().to_path_buf());
        let mut header = String::from("velocity_x,velocity_y");
        if *self.get_d() == 3 {
            header.push_str(",velocity_z");
        }
        if number_of_threads == 1 {
            let path = step_path.join(crate::io::VELOCITY_FILE);
            let mut file = File::create(path)?;
            writeln!(file, "{header}")?;
            self.get_nodes().iter().for_each(|node| {
                let line = node
                    .get_velocity()
                    .iter()
                    .map(|u_x| format!("{u_x:.8e}"))
                    .collect::<Vec<String>>()
                    .join(",");
                writeln!(file, "{line}").unwrap();
            });
        } else if number_of_threads > 1 {
            let number_of_nodes = self.get_number_of_nodes();
            let chunk_size = number_of_nodes / number_of_threads;
            self.get_nodes()
                .par_chunks(chunk_size)
                .enumerate()
                .for_each(|(i, chunk)| {
                    let file_str = format!("velocity_{i}.csv");
                    let path = step_path.join(file_str);
                    let mut file = File::create(path).unwrap();
                    writeln!(file, "{header}").unwrap();
                    chunk.iter().for_each(|node| {
                        let line = node
                            .get_velocity()
                            .iter()
                            .map(|u_x| format!("{u_x:.8e}"))
                            .collect::<Vec<String>>()
                            .join(",");
                        writeln!(file, "{line}").unwrap();
                    });
                });
        }
        Ok(())
    }

    pub fn write_coordinates(&self) -> LbResult<()> {
        let data_path = Path::new(crate::io::DATA_PATH);
        let path = data_path.join(crate::io::COORDINATES_FILE);
        println!("Writing {}.\n", crate::io::COORDINATES_FILE.yellow().bold());
        let mut file = File::create(path)?;
        match self.get_d() {
            2 => writeln!(file, "i,j,x,y,node_type")?,
            3 => writeln!(file, "i,j,k,x,y,z,node_type")?,
            _ => panic!("Unsupported dimension: {}", self.get_d()),
        }
        self.get_nodes().iter().for_each(|node| {
            let index = node
                .get_index()
                .iter()
                .map(|idx| idx.to_string())
                .collect::<Vec<String>>()
                .join(",");
            let coordinates = node
                .get_coordinates()
                .iter()
                .map(|coordinate| format!("{coordinate:.8e}"))
                .collect::<Vec<String>>()
                .join(",");
            let node_type = match node.get_node_type() {
                NodeType::Fluid => "fluid",
                NodeType::Solid => "solid",
            };
            writeln!(file, "{index},{coordinates},{node_type}").unwrap();
        });
        Ok(())
    }

    pub fn write_post_processing(&self, post_function: &super::post::PostFunction) -> LbResult<()> {
        if self.get_time_step() % post_function.interval == 0 {
            let post_results = &(post_function.function)(self);
            let post_processing_path = Path::new(crate::io::POST_PROCESSING_PATH);
            let path = post_processing_path.join(&post_function.file_name);
            let mut file = OpenOptions::new().create(true).append(true).open(path)?;
            if self.get_time_step() == 0 {
                write!(file, "time_step")?;
                for post_result in post_results {
                    write!(file, ",{}", post_result.name)?;
                }
                writeln!(file)?;
            }
            write!(file, "{}", self.get_time_step())?;
            for post_result in post_results {
                write!(file, ",{:.8e}", post_result.value)?;
            }
            writeln!(file)?;
        }
        Ok(())
    }
}

fn create_script_for_residuals_graph(dim: usize) -> LbResult<()> {
    let case_name = crate::io::get_case_name();
    let post_processing_path = Path::new(crate::io::POST_PROCESSING_PATH);
    let path = post_processing_path.join(crate::io::RESIDUALS_GRAPH_FILE);
    let mut file = File::create(&path)?;
    let path_str = path.to_str().unwrap();
    println!(
        "Creating the residuals graph script file: {}.\n",
        path_str.yellow().bold()
    );
    writeln!(
        file,
        r#"set datafile separator comma
set title "{case_name}"
set ylabel "Residuals"
set xlabel "Iterations"
set grid
set logscale y
set yrange [{min_tolerance}:]
set ytics format "%L"
set mxtics 5
set terminal push
set terminal pngcairo font "courier"
set output "fig_{case_name_prefix}_residuals.png"
plot "../data/residuals.csv" u 1:2 t "density" w l,\
"" u 1:3 t "velocity (x)" w l,\
"" u 1:4 t "velocity (y)" w l{velocity_z_str}
set terminal pdfcairo font "courier"
set output "fig_{case_name_prefix}_residuals.pdf"
replot
set terminal pop
set output"#,
        case_name = case_name.replace("_", "\\\\_"),
        min_tolerance = 1e-7,
        case_name_prefix = case_name.replace(" ", "_").to_lowercase(),
        velocity_z_str = if dim == 3 {
            ",\\\n\"\" u 1:5 t \"velocity (z)\" w l"
        } else {
            ""
        }
    )?;
    Ok(())
}

fn create_script_for_live_residuals_graph(dim: usize) -> LbResult<()> {
    let case_name = crate::io::get_case_name();
    let post_processing_path = Path::new(crate::io::POST_PROCESSING_PATH);
    let path = post_processing_path.join(crate::io::LIVE_RESIDUALS_GRAPH_FILE);
    let mut file = File::create(&path)?;
    let path_str = path.to_str().unwrap();
    println!(
        "Creating the live residuals graph script: {}.\n",
        path_str.yellow().bold()
    );
    writeln!(
        file,
        r#"set datafile separator comma
set title "{case_name}"
set ylabel "Residuals"
set xlabel "Iterations"
set grid
set logscale y
set yrange [{min_tolerance}:]
set ytics format "%L"
set mxtics 5
set terminal push
set terminal qt font "courier,12"
bind "q" "true=0"
true = 1
while (true) {{
plot "../data/residuals.csv" u 1:2 t "density" w l lw 2,\
"" u 1:3 t "velocity (x)" w l lw 2,\
"" u 1:4 t "velocity (y)" w l lw 2{velocity_z_str}
pause 1
}}
set terminal pop"#,
        case_name = case_name.replace("_", "\\\\_"),
        min_tolerance = 1e-7,
        velocity_z_str = if dim == 3 {
            ",\\\n\"\" u 1:5 t \"velocity (z)\" w l lw 2"
        } else {
            ""
        }
    )?;
    Ok(())
}

pub fn case_setup(momentum_parameters: &super::Parameters) {
    crate::io::create_case_directories().unwrap_or_else(|e| {
        eprintln! {"Error while creating the case directories: {e}"};
        std::process::exit(1);
    });
    let dim = momentum_parameters
        .velocity_set
        .get_velocity_set_parameters()
        .d;
    create_script_for_residuals_graph(dim).unwrap_or_else(|e| {
        eprintln! {"Error while creating the gnuplot script: {e}"};
        std::process::exit(1);
    });
    create_script_for_live_residuals_graph(dim).unwrap_or_else(|e| {
        eprintln! {"Error while creating the gnuplot script: {e}"};
        std::process::exit(1);
    });
}
