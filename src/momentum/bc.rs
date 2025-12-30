use super::Lattice;
use super::Node;
use crate::prelude_crate::*;
use rayon::prelude::*;

pub use BoundaryCondition::*;

pub(crate) const ZOUHE_ERROR_MESSAGE: &str =
    "Error: Wrong input of Zou & He (1997) boundary conditions.";

#[derive(Debug, PartialEq)]
pub enum BoundaryCondition {
    NoSlip,
    BounceBack {
        density: Float,
        velocity: Vec<Float>,
    },
    AntiBounceBack {
        density: Float,
    },
    ZouHe {
        density: Option<Float>,
        velocity: Vec<Option<Float>>,
    },
    Periodic,
}

impl Lattice {
    pub(super) fn boundary_conditions_step(&self) {
        self.get_boundary_nodes()
            .iter()
            .for_each(|(boundary_face, nodes)| {
                let boundary_condition = self.get_boundary_condition(boundary_face);
                match boundary_condition {
                    NoSlip => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_no_slip_bc(boundary_face);
                        });
                    }
                    BounceBack { density, velocity } => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_bounce_back_bc(boundary_face, density, velocity);
                        });
                    }
                    AntiBounceBack { density } => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_anti_bounce_back_bc(boundary_face, density);
                        });
                    }
                    Periodic => {}
                    ZouHe { density, velocity } => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_zou_he_bc(boundary_face, density, velocity);
                        });
                    }
                }
            });
    }
}

impl Node {
    fn compute_bounce_back_bc(
        &self,
        boundary_face: &BoundaryFace,
        density: &Float,
        velocity: &[Float],
    ) {
        let mut f = self.get_f();
        let f_star = self.get_f_star();
        let vel_set_params = self.get_vel_set_params();
        let w = vel_set_params.get_w();
        let c = vel_set_params.get_c();
        let q_faces = vel_set_params.get_q_faces(boundary_face);
        q_faces.iter().for_each(|&i| {
            let i_bar = vel_set_params.get_opposite_direction(i);
            let u_dot_c = velocity
                .iter()
                .zip(c[i].iter())
                .map(|(u_x, c_x)| u_x * (*c_x as Float))
                .sum::<Float>();
            f[i_bar] = f_star[i] - 2.0 * w[i] * density * CS_2_INV * u_dot_c;
        });
        self.set_f(f);
    }

    fn compute_anti_bounce_back_bc(&self, boundary_face: &BoundaryFace, density: &Float) {
        let mut f = self.get_f();
        let f_star = self.get_f_star();
        let vel_set_params = self.get_vel_set_params();
        let w = vel_set_params.get_w();
        let c = vel_set_params.get_c();
        let q_faces = vel_set_params.get_q_faces(boundary_face);
        let normal_direction = vel_set_params.get_face_normal_direction(boundary_face);
        let i_normal = normal_direction;
        let node_velocity = self.get_velocity();
        let neighbor_velocity = self.get_neighbor_node(i_normal).get_velocity();
        let velocity = node_velocity
            .iter()
            .zip(neighbor_velocity.iter())
            .map(|(u_x, nu_x)| 1.5 * u_x - 0.5 * nu_x)
            .collect::<Vec<Float>>();
        let u_dot_u = velocity.iter().map(|u_x| u_x * u_x).sum::<Float>();
        q_faces.iter().for_each(|&i| {
            let i_bar = vel_set_params.get_opposite_direction(i);
            let u_dot_c = velocity
                .iter()
                .zip(c[i].iter())
                .map(|(u_x, c_x)| u_x * (*c_x as Float))
                .sum::<Float>();
            f[i_bar] = -f_star[i]
                + 2.0
                    * w[i]
                    * density
                    * (1.0 + 0.5 * u_dot_c * u_dot_c * CS_4_INV - 0.5 * u_dot_u * CS_2_INV);
        });
        self.set_f(f);
    }

    fn compute_no_slip_bc(&self, boundary_face: &BoundaryFace) {
        let mut f = self.get_f();
        let f_star = self.get_f_star();
        let vel_set_params = self.get_vel_set_params();
        let q_faces = vel_set_params.get_q_faces(boundary_face);
        q_faces.iter().for_each(|&i| {
            let i_bar = vel_set_params.get_opposite_direction(i);
            f[i_bar] = f_star[i];
        });
        self.set_f(f);
    }

    fn compute_zou_he_bc(
        &self,
        boundary_face: &BoundaryFace,
        density: &Option<Float>,
        velocity: &[Option<Float>],
    ) {
        let f = self.get_f();
        let vel_set_params = self.get_vel_set_params();
        let f = vel_set_params
            .zou_he_bc_computation
            .as_ref()
            .expect("Zou & He BC computation function not defined.")(
            boundary_face,
            &f,
            density,
            velocity,
        );
        self.set_f(f);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_no_slip_bc_d2q9_west_face() {
        let node = Node::new(
            1.0,
            vec![0.0, 0.0],
            Arc::new(None),
            NodeType::Fluid,
            vec![0, 5],
            vec![0.005, 0.055],
            Arc::new(momentum::node::Parameters::test_default(2)),
        );
        node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
        node.set_f_star(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);

        node.compute_no_slip_bc(&BoundaryFace::West);

        assert_eq!(
            node.get_f(),
            vec![0.1, 4.0, 0.3, 0.4, 0.5, 8.0, 0.7, 0.8, 7.0]
        );
    }

    #[test]
    fn test_compute_no_slip_bc_d2q9_east_face() {
        let node = Node::new(
            1.0,
            vec![0.0, 0.0],
            Arc::new(None),
            NodeType::Fluid,
            vec![9, 5],
            vec![0.095, 0.055],
            Arc::new(momentum::node::Parameters::test_default(2)),
        );
        node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
        node.set_f_star(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);

        node.compute_no_slip_bc(&BoundaryFace::East);

        assert_eq!(
            node.get_f(),
            vec![0.1, 0.2, 0.3, 2.0, 0.5, 0.6, 9.0, 6.0, 0.9]
        );
    }

    #[test]
    fn test_compute_no_slip_bc_d2q9_south_face() {
        let node = Node::new(
            1.0,
            vec![0.0, 0.0],
            Arc::new(None),
            NodeType::Fluid,
            vec![9, 5],
            vec![0.095, 0.055],
            Arc::new(momentum::node::Parameters::test_default(2)),
        );
        node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
        node.set_f_star(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);

        node.compute_no_slip_bc(&BoundaryFace::South);

        assert_eq!(
            node.get_f(),
            vec![0.1, 0.2, 5.0, 0.4, 0.5, 8.0, 9.0, 0.8, 0.9]
        );
    }

    #[test]
    fn test_compute_no_slip_bc_d2q9_north_face() {
        let node = Node::new(
            1.0,
            vec![0.0, 0.0],
            Arc::new(None),
            NodeType::Fluid,
            vec![9, 5],
            vec![0.095, 0.055],
            Arc::new(momentum::node::Parameters::test_default(2)),
        );
        node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
        node.set_f_star(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);

        node.compute_no_slip_bc(&BoundaryFace::North);

        assert_eq!(
            node.get_f(),
            vec![0.1, 0.2, 0.3, 0.4, 3.0, 0.6, 0.7, 6.0, 7.0]
        );
    }
}
