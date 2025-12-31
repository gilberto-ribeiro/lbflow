use super::Lattice;
use super::Node;
use crate::prelude_crate::*;
use rayon::prelude::*;

pub use BoundaryCondition::*;
pub use InnerBoundaryCondition::*;

#[derive(Debug, PartialEq)]
pub enum BoundaryCondition {
    AntiBounceBack { scalar_value: Float },
    AntiBBNoFlux,
    BBNoFlux,
    ZerothOrderNoFlux,
    SecondOrderNoFlux,
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
                    ZerothOrderNoFlux => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_zeroth_order_no_flux_bc(boundary_face);
                        });
                    }
                    SecondOrderNoFlux => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_second_order_no_flux_bc(boundary_face);
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

    fn compute_zeroth_order_no_flux_bc(&self, boundary_face: &BoundaryFace) {
        // let mut g = self.get_g();
        let vel_set_params = self.get_vel_set_params();
        // let q_faces = vel_set_params.get_q_faces(boundary_face);
        let i_normal = vel_set_params.get_face_normal_direction(boundary_face);
        let neighbor_node = self.get_neighbor_node(i_normal);
        let neighbor_node_g = neighbor_node.get_g();
        // q_faces
        //     .iter()
        //     .map(|&i| vel_set_params.get_opposite_direction(i))
        //     .for_each(|i| g[i] = neighbor_g[i]);
        self.set_g(neighbor_node_g);
    }

    fn compute_second_order_no_flux_bc(&self, boundary_face: &BoundaryFace) {
        let mut g = self.get_g();
        let vel_set_params = self.get_vel_set_params();
        let q_faces = vel_set_params.get_q_faces(boundary_face);
        let i_normal = vel_set_params.get_face_normal_direction(boundary_face);
        let neighbor_node = self.get_neighbor_node(i_normal);
        let next_neighbor_node = neighbor_node.get_neighbor_node(i_normal);
        let neighbor_g = neighbor_node.get_g();
        let next_neighbor_g = next_neighbor_node.get_g();
        q_faces
            .iter()
            .map(|&i| vel_set_params.get_opposite_direction(i))
            .for_each(|i| g[i] = 2.0 * neighbor_g[i] - next_neighbor_g[i]);
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
