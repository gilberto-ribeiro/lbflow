use super::Lattice;
use crate::prelude::*;
use colored::*;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

impl Lattice {
    pub fn get_residuals_info(&self) -> crate::io::ResidualsInfo {
        let info = self.get_momentum_lattice().get_residuals_info();
        crate::io::ResidualsInfo {
            print_header: format!(
                "{} {:>16}",
                info.print_header,
                self.get_scalar_name().cyan().bold()
            ),
            print_line: format!(
                "{} {:>16.8e}",
                info.print_line,
                self.get_residuals().get_concentration()
            ),
            write_header: format!("{},{}", info.write_header, self.get_scalar_name()),
            write_line: format!(
                "{},{:.8e}",
                info.write_line,
                self.get_residuals().get_concentration()
            ),
            time_step: info.time_step,
        }
    }

    pub fn write_data(&self) {
        match self
            .get_momentum_lattice()
            .get_config()
            .get_write_data_mode()
        {
            WriteDataMode::Frequency(n) => {
                if self.get_momentum_lattice().get_time_step() % n == 0
                    || self.get_momentum_lattice().get_time_step() == 0
                {
                    println!();
                    self.write_data_from_steps();
                }
            }
            WriteDataMode::ListOfSteps(list) => {
                if list.contains(&self.get_momentum_lattice().get_time_step())
                    || self.get_momentum_lattice().get_time_step() == 0
                {
                    println!();
                    self.write_data_from_steps();
                }
            }
        }
    }

    fn write_data_from_steps(&self) {
        let data_path = Path::new(crate::io::DATA_PATH);
        let step_path = data_path.join(self.get_momentum_lattice().get_time_step().to_string());
        if let Err(e) = fs::create_dir_all(&step_path) {
            eprintln!(
                "Error while creating the the step {} directory: {e}.",
                self.get_momentum_lattice().get_time_step()
            );
            std::process::exit(1);
        };
        println!(
            "Writing {} for time step {}.\n",
            (self.get_scalar_name().to_string() + ".csv").yellow().bold(),
            self.get_momentum_lattice()
                .get_time_step()
                .to_string()
                .yellow()
                .bold()
        );
        if let Err(e) = self.write_concentration(&step_path) {
            eprintln!("Error while writing the concentration file: {e}.");
            std::process::exit(1);
        };
    }

    fn write_concentration<P>(&self, step_path: P) -> LbResult<()>
    where
        P: AsRef<Path>,
    {
        let number_of_threads = self
            .get_momentum_lattice()
            .get_config()
            .get_number_of_threads();
        let step_path = Arc::new(step_path.as_ref().to_path_buf());
        if number_of_threads == 1 {
            let path = step_path.join(self.get_scalar_name().to_string() + ".csv");
            let mut file = File::create(path)?;
            writeln!(file, "{}", self.get_scalar_name())?;
            self.get_nodes().iter().for_each(|node| {
                let line = format!("{:.8e}", node.get_concentration());
                writeln!(file, "{line}").unwrap();
            });
        } else if number_of_threads > 1 {
            let number_of_nodes = self.get_momentum_lattice().get_number_of_nodes();
            let chunk_size = number_of_nodes / number_of_threads;
            self.get_nodes()
                .par_chunks(chunk_size)
                .enumerate()
                .for_each(|(i, chunk)| {
                    let file_str = format!("{}_{}.csv", self.get_scalar_name(), i);
                    let path = step_path.join(file_str);
                    let mut file = File::create(path).unwrap();
                    writeln!(file, "{}", self.get_scalar_name()).unwrap();
                    chunk.iter().for_each(|node| {
                        let line = format!("{:.8e}", node.get_concentration());
                        writeln!(file, "{line}").unwrap();
                    });
                });
        }
        Ok(())
    }
}
