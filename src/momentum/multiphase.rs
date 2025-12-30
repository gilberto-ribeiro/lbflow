use super::Node;
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

impl Node {
    fn get_phi(&self) -> Option<Float> {
        *self.phi.read().unwrap()
    }

    fn set_phi(&self, phi: Float) {
        let mut phi_guard = self.phi.write().unwrap();
        *phi_guard = Some(phi);
    }

    fn get_multiphase_parameters(&self) -> &Option<multiphase::Parameters> {
        &self.multiphase_parameters
    }

    pub(super) fn compute_phi(&self) {
        if self.get_multiphase_parameters().is_some() {
            let density = self.get_density();
            let phi = kernel::pseudo_potential(density);
            self.set_phi(phi);
        }
    }

    pub(super) fn compute_shan_chen_force(&self) -> Option<Vec<Float>> {
        if let Some(mp_params) = self.get_multiphase_parameters() {
            let g_function = mp_params.get_g_function();
            let wall_phi = mp_params.get_wall_phi();
            let node_phi = self.get_phi().unwrap();
            let vs_params = self.get_vel_set_params();
            let c = vs_params.get_c();
            let w = vs_params.get_w();
            let d = vs_params.get_d();
            let q = vs_params.get_q();
            let mut neighbor_nodes_phi = vec![0.0; q];
            self.get_neighbor_nodes()
                .iter()
                .for_each(|(&i, neighbor_node)| {
                    neighbor_nodes_phi[i] = match neighbor_node.get_node_type() {
                        Fluid => neighbor_node.get_phi().unwrap(),
                        Solid => wall_phi,
                    };
                });
            let force = (0..d)
                .map(|x| {
                    -g_function
                        * node_phi
                        * (0..q)
                            .map(|i| {
                                let w_i = w[i];
                                let phi_i = neighbor_nodes_phi[i];
                                let c_ix = c[i][x] as Float;
                                w_i * c_ix * phi_i
                            })
                            .sum::<Float>()
                })
                .collect::<Vec<Float>>();
            Some(force)
        } else {
            None
        }
    }
}
