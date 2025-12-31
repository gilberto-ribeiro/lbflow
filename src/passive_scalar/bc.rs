use super::Lattice;
use super::Node;
use crate::prelude_crate::*;
use rayon::prelude::*;

pub use BoundaryCondition::*;
pub use InnerBoundaryCondition::*;

#[derive(Debug, PartialEq)]
pub enum BoundaryCondition {
    AntiBounceBack {
        scalar_value: Float,
    },
    AntiBBNoFlux,
    BBNoFlux,
    /// The option used in the papers below is `unknown_populations_only = false`
    /// - Klass, F., Gabbana, A., & Bartel, A. (2025). Perfectly Matched Layers and Characteristic Boundaries in Lattice Boltzmann: Accuracy vs Cost. **AIAA Journal**, 63(4), 1319–1329. https://doi.org/10.2514/1.J064563
    /// - Neeraj, T., Velten, C., Janiga, G., Zähringer, K., Namdar, R., Varnik, F., Thévenin, D., & Hosseini, S. A. (2023). Modeling Gas Flows in Packed Beds with the Lattice Boltzmann Method: Validation Against Experiments. **Flow, Turbulence and Combustion**, 111(2), 463–491. https://doi.org/10.1007/s10494-023-00444-z
    ZerothOrderNoFlux {
        unknown_populations_only: bool,
    },
    /// The option used in the papers below is `unknown_populations_only = true`
    /// - Junk, M., & Yang, Z. (2008). Outflow boundary conditions for the lattice Boltzmann method. **Progress in Computational Fluid Dynamics, An International Journal**, 8(1/2/3/4), 38. https://doi.org/10.1504/PCFD.2008.018077
    /// - Yu, D., Mei, R., & Shyy, W. (2005). Improved treatment of the open boundary in the method of Lattice Boltzmann equation: general description of the method. **Progress in Computational Fluid Dynamics, An International Journal**, 5(1/2), 3. https://doi.org/10.1504/PCFD.2005.005812
    SecondOrderNoFlux {
        unknown_populations_only: bool,
    },
    Periodic,
}

#[derive(Debug, PartialEq)]
pub enum InnerBoundaryCondition {
    InnerAntiBounceBack = 0,
    InnerBounceBack = 1,
}

impl<'a> Lattice<'a> {
    pub(super) fn boundary_conditions_step(&self) {
        self.get_boundary_nodes()
            .iter()
            .for_each(|(boundary_face, nodes)| {
                let boundary_condition = self.get_boundary_condition(boundary_face);
                let momentum_boundary_condition = self
                    .get_momentum_lattice()
                    .get_boundary_condition(boundary_face);
                let dim = *self.get_d();
                let velocity = match momentum_boundary_condition {
                    momentum::bc::NoSlip => Some(vec![0.0; dim]),
                    momentum::bc::BounceBack { velocity, .. } => Some(velocity.to_vec()),
                    momentum::bc::AntiBounceBack { .. } => None,
                    momentum::bc::Periodic => None,
                    momentum::bc::ZouHe { velocity, .. } => {
                        if velocity.iter().all(|v| v.is_some()) {
                            Some(velocity.iter().map(|v| v.unwrap()).collect::<Vec<Float>>())
                        } else {
                            None
                        }
                    }
                };
                match boundary_condition {
                    AntiBounceBack { scalar_value } => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_anti_bounce_back_bc(
                                boundary_face,
                                scalar_value,
                                velocity.as_deref(),
                            );
                        });
                    }
                    AntiBBNoFlux => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_anti_bb_no_flux_bc(boundary_face, velocity.as_deref());
                        });
                    }
                    BBNoFlux => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_bb_no_flux_bc(boundary_face);
                        });
                    }
                    ZerothOrderNoFlux {
                        unknown_populations_only,
                    } => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_zeroth_order_no_flux_bc(
                                boundary_face,
                                unknown_populations_only,
                            );
                        });
                    }
                    SecondOrderNoFlux {
                        unknown_populations_only,
                    } => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_second_order_no_flux_bc(
                                boundary_face,
                                unknown_populations_only,
                            );
                        });
                    }
                    Periodic => {}
                }
            });
    }
}

impl Node {
    fn compute_anti_bounce_back_bc(
        &self,
        boundary_face: &BoundaryFace,
        scalar_value: &Float,
        velocity: Option<&[Float]>,
    ) {
        let mut g = self.get_g();
        let g_star = self.get_g_star();
        let vel_set_params = self.get_vel_set_params();
        let w = vel_set_params.get_w();
        let c = vel_set_params.get_c();
        let q_faces = vel_set_params.get_q_faces(boundary_face);
        let velocity = self.get_wall_velocity(boundary_face, velocity);
        let u_dot_u = velocity.iter().map(|u_x| u_x * u_x).sum::<Float>();
        q_faces.iter().for_each(|&i| {
            let i_bar = vel_set_params.get_opposite_direction(i);
            let u_dot_c = velocity
                .iter()
                .zip(c[i].iter())
                .map(|(u_x, c_x)| u_x * (*c_x as Float))
                .sum::<Float>();
            let g_eq = w[i]
                * scalar_value
                * (1.0 + u_dot_c * CS_2_INV + 0.5 * u_dot_c * u_dot_c * CS_4_INV
                    - 0.5 * u_dot_u * CS_2_INV);
            g[i_bar] = -g_star[i] + 2.0 * g_eq;
        });
        self.set_g(g);
    }

    fn compute_anti_bb_no_flux_bc(&self, boundary_face: &BoundaryFace, velocity: Option<&[Float]>) {
        self.compute_anti_bounce_back_bc(boundary_face, &self.get_scalar_value(), velocity);
    }

    fn compute_bb_no_flux_bc(&self, boundary_face: &BoundaryFace) {
        let mut g = self.get_g();
        let g_star = self.get_g_star();
        let vel_set_params = self.get_vel_set_params();
        let q_faces = vel_set_params.get_q_faces(boundary_face);
        q_faces.iter().for_each(|&i| {
            let i_bar = vel_set_params.get_opposite_direction(i);
            g[i_bar] = g_star[i];
        });
        self.set_g(g);
    }

    fn compute_zeroth_order_no_flux_bc(
        &self,
        boundary_face: &BoundaryFace,
        unknown_populations_only: &bool,
    ) {
        let mut g = self.get_g();
        let vel_set_params = self.get_vel_set_params();
        let i_normal = vel_set_params.get_face_normal_direction(boundary_face);
        let neighbor_node = self.get_neighbor_node(i_normal);
        let neighbor_node_g = neighbor_node.get_g();
        if *unknown_populations_only {
            let q_faces = vel_set_params.get_q_faces(boundary_face);
            q_faces
                .iter()
                .map(|&i| vel_set_params.get_opposite_direction(i))
                .for_each(|i| g[i] = neighbor_node_g[i]);
        } else {
            g = neighbor_node_g;
        };
        self.set_g(g);
    }

    fn compute_second_order_no_flux_bc(
        &self,
        boundary_face: &BoundaryFace,
        unknown_populations_only: &bool,
    ) {
        let mut g = self.get_g();
        let vel_set_params = self.get_vel_set_params();
        let i_normal = vel_set_params.get_face_normal_direction(boundary_face);
        let neighbor_node = self.get_neighbor_node(i_normal);
        let next_neighbor_node = neighbor_node.get_neighbor_node(i_normal);
        let neighbor_g = neighbor_node.get_g();
        let next_neighbor_g = next_neighbor_node.get_g();
        if *unknown_populations_only {
            let q_faces = vel_set_params.get_q_faces(boundary_face);
            q_faces
                .iter()
                .map(|&i| vel_set_params.get_opposite_direction(i))
                .for_each(|i| g[i] = 2.0 * neighbor_g[i] - next_neighbor_g[i]);
        } else {
            let q = vel_set_params.get_q();
            (0..q).for_each(|i| g[i] = 2.0 * neighbor_g[i] - next_neighbor_g[i])
        };
        self.set_g(g);
    }

    fn get_wall_velocity(
        &self,
        boundary_face: &BoundaryFace,
        velocity: Option<&[Float]>,
    ) -> Vec<Float> {
        match velocity {
            Some(velocity) => velocity.to_vec(),
            None => {
                let vel_set_params = self.get_vel_set_params();
                let i_normal = vel_set_params.get_face_normal_direction(boundary_face);
                let node_velocity = self.get_velocity();
                let neighbor_velocity = self.get_neighbor_node(i_normal).get_velocity();
                node_velocity
                    .iter()
                    .zip(neighbor_velocity.iter())
                    .map(|(u_x, nu_x)| 1.5 * u_x - 0.5 * nu_x)
                    .collect::<Vec<Float>>()
            }
        }
    }
}
