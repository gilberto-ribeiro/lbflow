use crate::io::WriteDataMode;
use crate::parameters;

pub struct Parameters {
    pub case_name: String,
    pub data_write_mode: WriteDataMode,
}

pub struct Simulation {
    case_name: String,
    write_data_mode: WriteDataMode,
}

impl Simulation {
    pub fn new(parameters: Parameters) -> Self {
        Simulation {
            case_name: parameters.case_name,
            write_data_mode: parameters.data_write_mode,
        }
    }
}

impl Simulation {
    pub fn get_case_name(&self) -> &str {
        &self.case_name
    }

    pub fn get_data_write_mode(&self) -> &WriteDataMode {
        &self.write_data_mode
    }
}
