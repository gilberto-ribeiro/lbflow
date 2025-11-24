use super::Node;
use crate::prelude_crate::*;

pub use BoundaryCondition::*;

#[derive(Debug, PartialEq)]
pub enum BoundaryCondition {
    AntiBounceBack { scalar_value: Float },
    AntiBBNoFlux,
    Periodic,
}

impl Node {
    pub(super) fn compute_anti_bounce_back_bc(
        &self,
        boundary_face: &BoundaryFace,
        scalar_value: &Float,
        velocity: Option<&[Float]>,
    ) {
        let mut g = self.get_g();
        let g_star = self.get_g_star();
        let vel_set_params = self.get_velocity_set_parameters();
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

    pub(super) fn compute_no_flux_bc(
        &self,
        boundary_face: &BoundaryFace,
        velocity: Option<&[Float]>,
    ) {
        self.compute_anti_bounce_back_bc(boundary_face, &self.get_scalar_value(), velocity);
    }

    fn get_wall_velocity(
        &self,
        boundary_face: &BoundaryFace,
        velocity: Option<&[Float]>,
    ) -> Vec<Float> {
        match velocity {
            Some(velocity) => velocity.to_vec(),
            None => {
                let vel_set_params = self.get_velocity_set_parameters();
                let i_normal = vel_set_params.get_face_normal_direction(boundary_face);
                let node_velocity = self.get_momentum_node().get_velocity();
                let neighbor_velocity = self
                    .get_neighbor_node(i_normal)
                    .get_momentum_node()
                    .get_velocity();
                node_velocity
                    .iter()
                    .zip(neighbor_velocity.iter())
                    .map(|(u_x, nu_x)| 1.5 * u_x - 0.5 * nu_x)
                    .collect::<Vec<Float>>()
            }
        }
    }
}
