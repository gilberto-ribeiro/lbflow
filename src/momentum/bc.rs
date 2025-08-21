use super::Node;
use crate::prelude::*;

pub use BoundaryCondition::*;

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
    Periodic,
}

impl Node {
    pub fn compute_bounce_back_bc(
        &self,
        boundary_face: &BoundaryFace,
        density: &Float,
        velocity: &[Float],
    ) {
        let mut f = self.get_f();
        let f_star = self.get_f_star();
        let w = self.get_w();
        let c = self.get_c();
        let q_faces = self.get_q_faces(boundary_face);
        q_faces.iter().for_each(|&i| {
            let i_bar = self.get_opposite_direction(i);
            let u_dot_c = velocity
                .iter()
                .zip(c[i].iter())
                .map(|(u_x, c_x)| u_x * (*c_x as Float))
                .sum::<Float>();
            f[i_bar] = f_star[i] - 2.0 * w[i] * density * CS_2_INV * u_dot_c;
        });
        self.set_f(f);
    }

    pub fn compute_anti_bounce_back_bc(&self, boundary_face: &BoundaryFace, density: &Float) {
        let mut f = self.get_f();
        let f_star = self.get_f_star();
        let w = self.get_w();
        let c = self.get_c();
        let q_faces = self.get_q_faces(boundary_face);
        let normal_direction = self.get_face_normal_direction(boundary_face);
        let i_normal = *normal_direction;
        let node_velocity = self.get_velocity();
        let neighbor_velocity = self.get_neighbor_node(i_normal).get_velocity();
        let velocity = node_velocity
            .iter()
            .zip(neighbor_velocity.iter())
            .map(|(u_x, nu_x)| 1.5 * u_x - 0.5 * nu_x)
            .collect::<Vec<Float>>();
        let u_dot_u = velocity.iter().map(|u_x| u_x * u_x).sum::<Float>();
        q_faces.iter().for_each(|&i| {
            let i_bar = self.get_opposite_direction(i);
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

    /// # Example
    /// ## 2D
    /// ### West face
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// # use lbflow::momentum::ConversionFactor;
    /// # use lbflow::BoundaryFace;
    /// let node = Node::new(
    ///     1.0,
    ///     vec![0.0, 0.0],
    ///     NodeType::Fluid,
    ///     vec![0, 5],
    ///     vec![0.005, 0.055],
    ///     Arc::new(VelocitySetParameters::test_default(2)),
    ///     Arc::new(ConversionFactor::default()),
    /// );
    ///
    /// //                0    1    2    3    4    5    6    7    8
    /// node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    ///
    /// //                     0    1    2    3    4    5    6    7    8
    /// node.set_f_star(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    ///
    /// node.compute_no_slip_bc(&BoundaryFace::West);
    /// assert_eq!(node.get_f(), vec![0.1, 4.0, 0.3, 0.4, 0.5, 8.0, 0.7, 0.8, 7.0]);
    /// ```
    /// ### East face
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// # use lbflow::momentum::ConversionFactor;
    /// # use lbflow::BoundaryFace;
    /// let node = Node::new(
    ///     1.0,
    ///     vec![0.0, 0.0],
    ///     NodeType::Fluid,
    ///     vec![9, 5],
    ///     vec![0.095, 0.055],
    ///     Arc::new(VelocitySetParameters::test_default(2)),
    ///     Arc::new(ConversionFactor::default()),
    /// );
    ///
    /// //                0    1    2    3    4    5    6    7    8
    /// node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    ///
    /// //                     0    1    2    3    4    5    6    7    8
    /// node.set_f_star(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    ///
    /// node.compute_no_slip_bc(&BoundaryFace::East);
    /// assert_eq!(node.get_f(), vec![0.1, 0.2, 0.3, 2.0, 0.5, 0.6, 9.0, 6.0, 0.9]);
    /// ```
    /// ### South face
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// # use lbflow::momentum::ConversionFactor;
    /// # use lbflow::BoundaryFace;
    /// let node = Node::new(
    ///     1.0,
    ///     vec![0.0, 0.0],
    ///     NodeType::Fluid,
    ///     vec![9, 5],
    ///     vec![0.095, 0.055],
    ///     Arc::new(VelocitySetParameters::test_default(2)),
    ///     Arc::new(ConversionFactor::default()),
    /// );
    ///
    /// //                0    1    2    3    4    5    6    7    8
    /// node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    ///
    /// //                     0    1    2    3    4    5    6    7    8
    /// node.set_f_star(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    ///
    /// node.compute_no_slip_bc(&BoundaryFace::South);
    /// assert_eq!(node.get_f(), vec![0.1, 0.2, 5.0, 0.4, 0.5, 8.0, 9.0, 0.8, 0.9]);
    /// ```
    /// ### North face
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// # use lbflow::momentum::ConversionFactor;
    /// # use lbflow::BoundaryFace;
    /// let node = Node::new(
    ///     1.0,
    ///     vec![0.0, 0.0],
    ///     NodeType::Fluid,
    ///     vec![9, 5],
    ///     vec![0.095, 0.055],
    ///     Arc::new(VelocitySetParameters::test_default(2)),
    ///     Arc::new(ConversionFactor::default()),
    /// );
    ///
    /// //                0    1    2    3    4    5    6    7    8
    /// node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    ///
    /// //                     0    1    2    3    4    5    6    7    8
    /// node.set_f_star(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    ///
    /// node.compute_no_slip_bc(&BoundaryFace::North);
    /// assert_eq!(node.get_f(), vec![0.1, 0.2, 0.3, 0.4, 3.0, 0.6, 0.7, 6.0, 7.0]);
    /// ```
    pub fn compute_no_slip_bc(&self, boundary_face: &BoundaryFace) {
        let mut f = self.get_f();
        let f_star = self.get_f_star();
        let q_faces = self.get_q_faces(boundary_face);
        q_faces.iter().for_each(|&i| {
            let i_bar = self.get_opposite_direction(i);
            f[i_bar] = f_star[i];
        });
        self.set_f(f);
    }
}
