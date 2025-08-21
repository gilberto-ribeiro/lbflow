use super::Node;
use crate::prelude::*;

pub use BoundaryCondition::*;

#[derive(Debug, PartialEq)]
pub enum BoundaryCondition {
    AntiBounceBack { concentration: Float },
    AntiBBNoFlux,
    Periodic,
}

impl Node {
    pub fn compute_anti_bounce_back_bc(
        &self,
        boundary_face: &BoundaryFace,
        concentration: &Float,
        velocity: Option<&[Float]>,
    ) {
        let mut g = self.get_g();
        let g_star = self.get_g_star();
        let w = self.get_w();
        let c = self.get_c();
        let q_faces = self.get_q_faces(boundary_face);
        let velocity = self.get_wall_velocity(boundary_face, velocity);
        let u_dot_u = velocity.iter().map(|u_x| u_x * u_x).sum::<Float>();
        q_faces.iter().for_each(|&i| {
            let i_bar = self.get_opposite_direction(i);
            let u_dot_c = velocity
                .iter()
                .zip(c[i].iter())
                .map(|(u_x, c_x)| u_x * (*c_x as Float))
                .sum::<Float>();
            let g_eq = w[i]
                * concentration
                * (1.0 + u_dot_c * CS_2_INV + 0.5 * u_dot_c * u_dot_c * CS_4_INV
                    - 0.5 * u_dot_u * CS_2_INV);
            g[i_bar] = -g_star[i] + 2.0 * g_eq;
        });
        self.set_g(g);
    }

    pub fn compute_no_flux_bc(&self, boundary_face: &BoundaryFace, velocity: Option<&[Float]>) {
        self.compute_anti_bounce_back_bc(boundary_face, &self.get_concentration(), velocity);
    }

    fn get_wall_velocity(
        &self,
        boundary_face: &BoundaryFace,
        velocity: Option<&[Float]>,
    ) -> Vec<Float> {
        match velocity {
            Some(velocity) => velocity.to_vec(),
            None => {
                let i_normal = *self.get_face_normal_direction(boundary_face);
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
