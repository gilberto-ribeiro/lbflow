use super::Lattice;
use crate::NodeType;
use crate::io::WriteDataMode;
use colored::*;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::Path;
use std::process;

impl Lattice {
    pub fn print_residuals(&self) {
        let residuals_density = self.get_residuals().get_density();
        let binding = self.get_residuals();
        let residuals_velocity = binding.get_velocity();
        if self.get_time_step() % 100 == 0 {
            let directions = ["x", "y", "z"];
            let mut header = format!(
                "\n{:>10} {:>16}",
                "time_step".cyan().bold(),
                "density".cyan().bold()
            );
            residuals_velocity
                .iter()
                .zip(directions.iter())
                .for_each(|(_, x)| {
                    header.push_str(&format!(" {:>16}", format!("velocity_{x}").cyan().bold()));
                });
            println!("{header}\n");
        }
        let mut line = format! {"{:>10} {:>16.8e}", self.get_time_step(), residuals_density};
        residuals_velocity.iter().for_each(|u_x| {
            line.push_str(&format!(" {u_x:>16.8e}"));
        });
        println!("{line}");
    }

    pub fn write_residuals(&self) -> io::Result<()> {
        let residuals_density = self.get_residuals().get_density();
        let binding = self.get_residuals();
        let residuals_velocity = binding.get_velocity();
        let data_path = Path::new(crate::io::DATA_PATH);
        let path = data_path.join(crate::io::RESIDUALS_FILE);
        let mut file = OpenOptions::new().create(true).append(true).open(path)?;
        if self.get_time_step() == 0 {
            let mut header = String::from("time_step,density");
            let directions = ["x", "y", "z"];
            let velocities = (0..*self.get_d())
                .map(|x| format!("velocity_{}", directions[x]))
                .collect::<Vec<String>>()
                .join(",");
            header = header + "," + &velocities;
            writeln!(file, "{header}")?;
        }
        let mut line = format! {"{},{:.8e}", self.get_time_step(), residuals_density};
        residuals_velocity.iter().for_each(|u_x| {
            line.push_str(&format!(",{u_x:.8e}"));
        });
        writeln!(file, "{line}")?;
        Ok(())
    }

    pub fn write_data(&self, write_data_mode: &WriteDataMode) {
        match write_data_mode {
            WriteDataMode::Frequency(n) => {
                if self.get_time_step() % n == 0 || self.get_time_step() == 0 {
                    println!();
                    self.write_data_from_steps();
                    self.write_vtk_from_steps();
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

    pub fn write_data_from_steps(&self) {
        let data_path = Path::new(crate::io::DATA_PATH);
        let step_path = data_path.join(self.get_time_step().to_string());
        if let Err(e) = fs::create_dir_all(&step_path) {
            eprintln!(
                "Error while creating the the step {} directory: {e}.",
                self.get_time_step()
            );
            process::exit(1);
        };
        let path = step_path.join(crate::io::DENSITY_FILE);
        println!(
            "\nWriting {} for time step {}.\n",
            crate::io::DENSITY_FILE.yellow().bold(),
            self.get_time_step().to_string().yellow().bold()
        );
        if let Err(e) = self.write_density(path) {
            eprintln!("Error while writing the density file: {e}.");
            process::exit(1);
        };
        let path = step_path.join(crate::io::VELOCITY_FILE);
        println!(
            "\nWriting {} for time step {}.\n",
            crate::io::VELOCITY_FILE.yellow().bold(),
            self.get_time_step().to_string().yellow().bold(),
        );
        if let Err(e) = self.write_velocity(path) {
            eprintln!("Error while writing the velocity file: {e}.");
            process::exit(1);
        };
    }

    fn write_density<P>(&self, path: P) -> io::Result<()>
    where
        P: AsRef<Path>,
    {
        let mut file = File::create(path)?;
        writeln!(file, "density")?;
        self.get_nodes().iter().for_each(|node| {
            let line = format!("{:.8e}", node.get_density());
            writeln!(file, "{line}").unwrap();
        });
        Ok(())
    }

    fn write_velocity<P>(&self, path: P) -> io::Result<()>
    where
        P: AsRef<Path>,
    {
        let mut file = File::create(path)?;
        let mut header = String::from("velocity_x,velocity_y");
        if *self.get_d() == 3 {
            header.push_str(",velocity_z");
        }
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
        Ok(())
    }

    pub fn write_coordinates(&self) -> io::Result<()> {
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
            let line = index + "," + &coordinates + "," + &node_type;
            writeln!(file, "{line}").unwrap();
        });
        Ok(())
    }

    pub fn write_vtk_from_steps(&self) {
        let vtk_files_path = Path::new(crate::io::VTK_PATH);
        let case_name = "case_test".replace(" ", "_").to_lowercase();
        let path_str = format!("{}_results_{:08}.vtk", case_name, self.get_time_step());
        let path = vtk_files_path.join(&path_str);
        println!(
            "\nWriting {} for time step {}.\n",
            path_str.yellow().bold(),
            self.get_time_step().to_string().yellow().bold()
        );
        if let Err(e) = self.write_vtk(path) {
            eprintln!("Error while writing the vtk file: {e}.");
            process::exit(1);
        };
    }

    fn write_vtk<P>(&self, path: P) -> io::Result<()>
    where
        P: AsRef<Path>,
    {
        let point_data = self.get_number_of_nodes();
        let mut file = File::create(path)?;
        writeln!(file, "# vtk DataFile Version 3.0")?;
        writeln!(file, "LBM simulation data")?;
        writeln!(file, "ASCII")?;
        writeln!(file, "DATASET STRUCTURED_GRID")?;
        writeln!(
            file,
            "DIMENSIONS {} {} {}",
            self.get_nx(),
            self.get_ny(),
            self.get_nz()
        )?;
        writeln!(file, "POINTS {point_data} float")?;
        self.get_nodes().iter().for_each(|node| {
            let line = node
                .get_coordinates()
                .iter()
                .map(|coord| format!("{coord:.4e}"))
                .collect::<Vec<String>>()
                .join(" ");
            writeln!(file, "{line}").unwrap();
        });
        writeln!(file, "POINT_DATA {point_data}")?;
        writeln!(file, "SCALARS density float 1")?;
        writeln!(file, "LOOKUP_TABLE default")?;
        self.get_nodes().iter().for_each(|node| {
            writeln!(file, "{:.8e}", node.get_density()).unwrap();
        });
        writeln!(file, "VECTORS velocity float")?;
        self.nodes.iter().for_each(|node| {
            let line = node
                .get_velocity()
                .iter()
                .map(|u_x| format!("{u_x:.8e}"))
                .collect::<Vec<String>>()
                .join(" ");
            writeln!(file, "{line}").unwrap();
        });
        Ok(())
    }
}
