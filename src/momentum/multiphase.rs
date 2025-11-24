use crate::prelude_crate::*;

pub struct Parameters {
    g_function: Float,
    wall_density: Float,
    wall_phi: Float,
}

impl Parameters {
    pub fn new(g_function: Float, wall_density: Float) -> Self {
        let wall_phi = kernel::pseudo_potential(wall_density);
        Parameters {
            g_function,
            wall_density,
            wall_phi,
        }
    }
}

impl Parameters {
    pub(super) fn get_g_function(&self) -> Float {
        self.g_function
    }

    fn _get_wall_density(&self) -> Float {
        self.wall_density
    }

    pub(super) fn get_wall_phi(&self) -> Float {
        self.wall_phi
    }
}
