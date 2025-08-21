use super::ConversionFactor;
use crate::prelude::*;
use crate::velocity_set::VelocityComputation;

// -------------------------------------------------------------------------- STRUCT: Node

#[derive(Debug)]
pub struct Node {
    density: RwLock<Float>,
    velocity: RwLock<Vec<Float>>,
    f: RwLock<Vec<Float>>,
    f_eq: RwLock<Vec<Float>>,
    f_star: RwLock<Vec<Float>>,
    node_type: NodeType,
    index: Vec<usize>,
    coordinates: Vec<Float>,
    velocity_set_parameters: Arc<VelocitySetParameters>,
    conversion_factor: Arc<ConversionFactor>,
    neighbor_nodes: RwLock<Option<HashMap<usize, Arc<Node>>>>,
    bounce_back_neighbor_nodes: RwLock<Option<HashMap<usize, Arc<Node>>>>,
    shallow_node: ShallowNode,
}

impl Node {
    pub fn new(
        density: Float,
        velocity: Vec<Float>,
        node_type: NodeType,
        index: Vec<usize>,
        coordinates: Vec<Float>,
        velocity_set_parameters: Arc<VelocitySetParameters>,
        conversion_factor: Arc<ConversionFactor>,
    ) -> Self {
        let q = velocity_set_parameters.q;
        Node {
            density: RwLock::new(density),
            velocity: RwLock::new(velocity.clone()),
            f: RwLock::new(vec![0.0; q]),
            f_eq: RwLock::new(vec![0.0; q]),
            f_star: RwLock::new(vec![0.0; q]),
            node_type,
            index,
            coordinates,
            velocity_set_parameters,
            conversion_factor,
            neighbor_nodes: RwLock::new(None),
            bounce_back_neighbor_nodes: RwLock::new(None),
            shallow_node: ShallowNode::new(density, velocity.clone()),
        }
    }

    pub fn test_default(dim: usize) -> Self {
        match dim {
            2 => Node::new(
                1.0,
                vec![0.0, 0.0],
                Fluid,
                vec![3, 7],
                vec![0.035, 0.075],
                Arc::new(VelocitySetParameters::test_default(2)),
                Arc::new(ConversionFactor::default()),
            ),
            3 => Node::new(
                1.0,
                vec![0.0, 0.0, 0.0],
                Fluid,
                vec![3, 5, 7],
                vec![0.035, 0.055, 0.075],
                Arc::new(VelocitySetParameters::test_default(3)),
                Arc::new(ConversionFactor::default()),
            ),
            _ => panic!("Unsupported dimension: {dim}"),
        }
    }
}

impl Default for Node {
    fn default() -> Self {
        Node::new(
            1.0,
            vec![0.0, 0.0],
            NodeType::Fluid,
            vec![0, 0],
            vec![0.0, 0.0],
            Arc::new(VelocitySetParameters::default()),
            Arc::new(ConversionFactor::default()),
        )
    }
}

impl Node {
    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_density(), 1.0);
    /// assert_eq!(&node.get_density(), &1.0);
    /// ```
    pub fn get_density(&self) -> Float {
        // *self.density.borrow()
        *self.density.read().unwrap()
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_density(), 1.0);
    ///
    /// node.set_density(1.1);
    /// assert_eq!(node.get_density(), 1.1);
    ///
    /// node.set_density(1.2);
    /// assert_eq!(node.get_density(), 1.2);
    /// ```
    fn set_density(&self, density: Float) {
        // self.density.replace(density);
        let mut density_guard = self.density.write().unwrap();
        *density_guard = density;
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_velocity(), vec![0.0, 0.0]);
    /// assert_eq!(node.get_velocity(), [0.0, 0.0]);
    /// assert_eq!(node.get_velocity(), &[0.0, 0.0]);
    /// assert_eq!(&node.get_velocity(), &[0.0, 0.0]);
    /// ```
    /// ## 3D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(3);
    ///
    /// assert_eq!(node.get_velocity(), vec![0.0, 0.0, 0.0]);
    /// assert_eq!(node.get_velocity(), [0.0, 0.0, 0.0]);
    /// assert_eq!(node.get_velocity(), &[0.0, 0.0, 0.0]);
    /// assert_eq!(&node.get_velocity(), &[0.0, 0.0, 0.0]);
    /// ```
    pub fn get_velocity(&self) -> Vec<Float> {
        // self.velocity.borrow().clone()
        self.velocity.read().unwrap().clone()
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_velocity(), vec![0.0, 0.0]);
    ///
    /// node.set_velocity(vec![0.1, 0.0]);
    /// let mut velocity = node.get_velocity();
    /// velocity[1] = 0.2;
    /// node.set_velocity(velocity);
    /// assert_eq!(node.get_velocity(), vec![0.1, 0.2]);
    ///
    /// let mut velocity = Vec::with_capacity(2);
    /// let mut u_x = 0.3;
    /// for i in 0..2 {
    ///     velocity.push(u_x);
    ///     u_x += 0.1;
    /// }
    /// node.set_velocity(velocity);
    /// let actual = node.get_velocity();
    /// let target = vec![0.3, 0.4];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    /// let mut velocity = vec![0.0; 2];
    /// let mut u_x = 0.5;
    /// for i in 0..2 {
    ///     velocity[i] = u_x;
    ///     u_x += 0.1;
    /// }
    /// node.set_velocity(velocity);
    /// let actual = node.get_velocity();
    /// let target = vec![0.5, 0.6];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    /// ```
    /// # 3D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(3);
    ///
    /// assert_eq!(node.get_velocity(), vec![0.0, 0.0, 0.0]);
    ///
    /// node.set_velocity(vec![0.1, 0.2, 0.0]);
    /// let mut velocity = node.get_velocity();
    /// velocity[2] = 0.3;
    /// node.set_velocity(velocity);
    /// assert_eq!(node.get_velocity(), vec![0.1, 0.2, 0.3]);
    /// ```
    fn set_velocity(&self, velocity: Vec<Float>) {
        // self.velocity.replace(velocity);
        let mut velocity_guard = self.velocity.write().unwrap();
        *velocity_guard = velocity;
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_f(), vec![0.0; 9]);
    /// assert_eq!(node.get_f(), [0.0; 9]);
    /// assert_eq!(node.get_f(), &[0.0; 9]);
    /// assert_eq!(&node.get_f(), &[0.0; 9]);
    /// ```
    pub fn get_f(&self) -> Vec<Float> {
        // self.f.borrow().clone()
        self.f.read().unwrap().clone()
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::constants::Float;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_f(), vec![0.0; 9]);
    ///
    /// node.set_f(vec![0.1, 0.0, 0.3, 0.0, 0.5, 0.0, 0.7, 0.0, 0.9]);
    /// let mut f = node.get_f();
    /// f[1] = 0.2;
    /// f[3] = 0.4;
    /// f[5] = 0.6;
    /// f[7] = 0.8;
    /// node.set_f(f);
    /// assert_eq!(node.get_f(), vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    ///
    /// let mut f = Vec::with_capacity(9);
    /// let mut f_i = 0.9;
    /// for i in 0..9 {
    ///    f.push(f_i);
    ///    f_i -= 0.1;
    /// }
    /// node.set_f(f);
    /// let actual = node.get_f();
    /// let target = vec![0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    /// let mut f = vec![0.0; 9];
    /// let mut f_i = 0.05;
    /// for i in 0..9 {
    ///    f[i] = f_i;
    ///    f_i += 0.05;
    /// }
    /// node.set_f(f);
    /// let actual = node.get_f();
    /// let target = vec![0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    ///
    /// let node = Node::test_default(3);
    ///
    /// assert_eq!(node.get_f(), vec![0.0; 27]);
    ///
    /// node.set_f(vec![0.1, 0.2, 0.0, 0.4, 0.5, 0.6, 0.7, 0.6, 0.0, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.0, 0.6, 0.7, 0.6, 0.5, 0.4, 0.0, 0.2, 0.1]);
    /// let mut f = node.get_f();
    /// f[2] = 0.3;
    /// f[8] = 0.5;
    /// f[18] = 0.5;
    /// f[24] = 0.3;
    /// node.set_f(f);
    /// assert_eq!(node.get_f(), vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]);
    /// ```
    pub fn set_f(&self, f: Vec<Float>) {
        // self.f.replace(f);
        let mut f_guard = self.f.write().unwrap();
        *f_guard = f;
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_f_eq(), vec![0.0; 9]);
    /// assert_eq!(node.get_f_eq(), [0.0; 9]);
    /// assert_eq!(node.get_f_eq(), &[0.0; 9]);
    /// assert_eq!(&node.get_f_eq(), &[0.0; 9]);
    /// ```
    pub fn get_f_eq(&self) -> Vec<Float> {
        // self.f_eq.borrow().clone()
        self.f_eq.read().unwrap().clone()
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::constants::Float;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_f_eq(), vec![0.0; 9]);
    ///
    /// node.set_f_eq(vec![0.1, 0.0, 0.3, 0.0, 0.5, 0.0, 0.7, 0.0, 0.9]);
    /// let mut f_eq = node.get_f_eq();
    /// f_eq[1] = 0.2;
    /// f_eq[3] = 0.4;
    /// f_eq[5] = 0.6;
    /// f_eq[7] = 0.8;
    /// node.set_f_eq(f_eq);
    /// assert_eq!(node.get_f_eq(), vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    ///
    /// let mut f_eq = Vec::with_capacity(9);
    /// let mut f_eq_i = 0.9;
    /// for i in 0..9 {
    ///    f_eq.push(f_eq_i);
    ///    f_eq_i -= 0.1;
    /// }
    /// node.set_f_eq(f_eq);
    /// let actual = node.get_f_eq();
    /// let target = vec![0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    /// let mut f_eq = vec![0.0; 9];
    /// let mut f_eq_i = 0.05;
    /// for i in 0..9 {
    ///    f_eq[i] = f_eq_i;
    ///    f_eq_i += 0.05;
    /// }
    /// node.set_f_eq(f_eq);
    /// let actual = node.get_f_eq();
    /// let target = vec![0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    ///
    /// let node = Node::test_default(3);
    ///
    /// assert_eq!(node.get_f_eq(), vec![0.0; 27]);
    ///
    /// node.set_f_eq(vec![0.1, 0.2, 0.0, 0.4, 0.5, 0.6, 0.7, 0.6, 0.0, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.0, 0.6, 0.7, 0.6, 0.5, 0.4, 0.0, 0.2, 0.1]);
    /// let mut f_eq = node.get_f_eq();
    /// f_eq[2] = 0.3;
    /// f_eq[8] = 0.5;
    /// f_eq[18] = 0.5;
    /// f_eq[24] = 0.3;
    /// node.set_f_eq(f_eq);
    /// assert_eq!(node.get_f_eq(), vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]);
    /// ```
    fn set_f_eq(&self, f_eq: Vec<Float>) {
        // self.f_eq.replace(f_eq);
        let mut f_eq_guard = self.f_eq.write().unwrap();
        *f_eq_guard = f_eq;
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_f_star(), vec![0.0; 9]);
    /// assert_eq!(node.get_f_star(), [0.0; 9]);
    /// assert_eq!(node.get_f_star(), &[0.0; 9]);
    /// assert_eq!(&node.get_f_star(), &[0.0; 9]);
    /// ```
    pub fn get_f_star(&self) -> Vec<Float> {
        // self.f_star.borrow().clone()
        self.f_star.read().unwrap().clone()
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::constants::Float;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_f_star(), vec![0.0; 9]);
    ///
    /// node.set_f_star(vec![0.1, 0.0, 0.3, 0.0, 0.5, 0.0, 0.7, 0.0, 0.9]);
    /// let mut f_star = node.get_f_star();
    /// f_star[1] = 0.2;
    /// f_star[3] = 0.4;
    /// f_star[5] = 0.6;
    /// f_star[7] = 0.8;
    /// node.set_f_star(f_star);
    /// assert_eq!(node.get_f_star(), vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    ///
    /// let mut f_star = Vec::with_capacity(9);
    /// let mut f_star_i = 0.9;
    /// for i in 0..9 {
    ///    f_star.push(f_star_i);
    ///    f_star_i -= 0.1;
    /// }
    /// node.set_f_star(f_star);
    /// let actual = node.get_f_star();
    /// let target = vec![0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    /// let mut f_star = vec![0.0; 9];
    /// let mut f_star_i = 0.05;
    /// for i in 0..9 {
    ///    f_star[i] = f_star_i;
    ///    f_star_i += 0.05;
    /// }
    /// node.set_f_star(f_star);
    /// let actual = node.get_f_star();
    /// let target = vec![0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    ///
    /// let node = Node::test_default(3);
    ///
    /// assert_eq!(node.get_f_star(), vec![0.0; 27]);
    ///
    /// node.set_f_star(vec![0.1, 0.2, 0.0, 0.4, 0.5, 0.6, 0.7, 0.6, 0.0, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.0, 0.6, 0.7, 0.6, 0.5, 0.4, 0.0, 0.2, 0.1]);
    /// let mut f_star = node.get_f_star();
    /// f_star[2] = 0.3;
    /// f_star[8] = 0.5;
    /// f_star[18] = 0.5;
    /// f_star[24] = 0.3;
    /// node.set_f_star(f_star);
    /// assert_eq!(node.get_f_star(), vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]);
    /// ```
    fn set_f_star(&self, f_star: Vec<Float>) {
        // self.f_star.replace(f_star);
        let mut f_star_guard = self.f_star.write().unwrap();
        *f_star_guard = f_star;
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_node_type(), &NodeType::Fluid);
    /// assert_ne!(node.get_node_type(), &NodeType::Solid);
    /// ```
    pub fn get_node_type(&self) -> &NodeType {
        &self.node_type
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_index(), &vec![3, 7]);
    /// assert_eq!(node.get_index(), &[3, 7]);
    /// ```
    pub fn get_index(&self) -> &Vec<usize> {
        &self.index
    }

    pub fn get_coordinates(&self) -> &Vec<Float> {
        &self.coordinates
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(2);
    ///
    /// let node = momentum_lattice.get_node_by_index(&vec![3, 7]);
    /// let neighbor_nodes = node.get_neighbor_nodes();
    ///
    /// assert_eq!(neighbor_nodes[&0].get_index(), &vec![3, 7]);
    /// assert_eq!(neighbor_nodes[&1].get_index(), &vec![4, 7]);
    /// assert_eq!(neighbor_nodes[&2].get_index(), &vec![3, 8]);
    /// assert_eq!(neighbor_nodes[&3].get_index(), &vec![2, 7]);
    /// assert_eq!(neighbor_nodes[&4].get_index(), &vec![3, 6]);
    /// assert_eq!(neighbor_nodes[&5].get_index(), &vec![4, 8]);
    /// assert_eq!(neighbor_nodes[&6].get_index(), &vec![2, 8]);
    /// assert_eq!(neighbor_nodes[&7].get_index(), &vec![2, 6]);
    /// assert_eq!(neighbor_nodes[&8].get_index(), &vec![4, 6]);
    /// ```
    /// ## 3D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(3);
    ///
    /// let node = momentum_lattice.get_node_by_index(&vec![3, 5, 7]);
    /// let neighbor_nodes = node.get_neighbor_nodes();
    ///
    /// assert_eq!(neighbor_nodes[&0].get_index(), &vec![3, 5, 7]);
    /// assert_eq!(neighbor_nodes[&1].get_index(), &vec![4, 5, 7]);
    /// assert_eq!(neighbor_nodes[&7].get_index(), &vec![4, 6, 7]);
    /// assert_eq!(neighbor_nodes[&19].get_index(), &vec![4, 6, 8]);
    /// ```
    pub fn get_neighbor_nodes(&self) -> HashMap<usize, Arc<Node>> {
        // self.neighbor_nodes.borrow().as_ref().cloned().unwrap()
        self.neighbor_nodes
            .read()
            .unwrap()
            .as_ref()
            .cloned()
            .unwrap()
    }

    pub(super) fn set_neighbor_nodes(&self, neighbor_nodes: HashMap<usize, Arc<Node>>) {
        // self.neighbor_nodes.replace(Some(neighbors));
        let mut neighbor_nodes_guard = self.neighbor_nodes.write().unwrap();
        *neighbor_nodes_guard = Some(neighbor_nodes);
    }

    pub fn get_bounce_back_neighbor_nodes(&self) -> HashMap<usize, Arc<Node>> {
        // self.bounce_back_neighbor_nodes
        //     .borrow()
        //     .as_ref()
        //     .cloned()
        //     .unwrap()
        self.bounce_back_neighbor_nodes
            .read()
            .unwrap()
            .as_ref()
            .cloned()
            .unwrap()
    }

    pub(super) fn set_bounce_back_neighbor_nodes(
        &self,
        bounce_back_neighbor_nodes: HashMap<usize, Arc<Node>>,
    ) {
        // node.bounce_back_neighbor_nodes
        //     .replace(Some(bounce_back_neighbor_nodes));
        let mut bounce_back_neighbor_nodes_guard = self.bounce_back_neighbor_nodes.write().unwrap();
        *bounce_back_neighbor_nodes_guard = Some(bounce_back_neighbor_nodes);
    }

    /// # Examples
    /// ## 2D
    /// ### Inner node
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(2);
    ///
    /// let node = momentum_lattice.get_node_by_index(&vec![3, 7]);
    ///
    /// assert_eq!(node.get_neighbor_node(0).get_index(), &vec![3, 7]);
    /// assert_eq!(node.get_neighbor_node(1).get_index(), &vec![4, 7]);
    /// assert_eq!(node.get_neighbor_node(2).get_index(), &vec![3, 8]);
    /// assert_eq!(node.get_neighbor_node(3).get_index(), &vec![2, 7]);
    /// assert_eq!(node.get_neighbor_node(4).get_index(), &vec![3, 6]);
    /// assert_eq!(node.get_neighbor_node(5).get_index(), &vec![4, 8]);
    /// assert_eq!(node.get_neighbor_node(6).get_index(), &vec![2, 8]);
    /// assert_eq!(node.get_neighbor_node(7).get_index(), &vec![2, 6]);
    /// assert_eq!(node.get_neighbor_node(8).get_index(), &vec![4, 6]);
    /// ```
    /// ### Edge node
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(2);
    ///
    /// let node = momentum_lattice.get_node_by_index(&vec![9, 9]);
    ///
    /// assert_eq!(node.get_neighbor_node(0).get_index(), &vec![9, 9]);
    /// assert_eq!(node.get_neighbor_node(1).get_index(), &vec![0, 9]);
    /// assert_eq!(node.get_neighbor_node(2).get_index(), &vec![9, 0]);
    /// assert_eq!(node.get_neighbor_node(3).get_index(), &vec![8, 9]);
    /// assert_eq!(node.get_neighbor_node(4).get_index(), &vec![9, 8]);
    /// assert_eq!(node.get_neighbor_node(5).get_index(), &vec![0, 0]);
    /// assert_eq!(node.get_neighbor_node(6).get_index(), &vec![8, 0]);
    /// assert_eq!(node.get_neighbor_node(7).get_index(), &vec![8, 8]);
    /// assert_eq!(node.get_neighbor_node(8).get_index(), &vec![0, 8]);
    ///
    /// assert_eq!(node.get_neighbor_node(3).get_neighbor_node(2).get_index(), &vec![8, 0]);
    /// ```
    /// ## 3D
    /// ### Inner node
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(3);
    ///
    /// let node = momentum_lattice.get_node_by_index(&vec![3, 5, 7]);
    ///
    /// assert_eq!(node.get_neighbor_node(0).get_index(), &vec![3, 5, 7]);
    /// assert_eq!(node.get_neighbor_node(1).get_index(), &vec![4, 5, 7]);
    /// assert_eq!(node.get_neighbor_node(7).get_index(), &vec![4, 6, 7]);
    /// assert_eq!(node.get_neighbor_node(19).get_index(), &vec![4, 6, 8]);
    /// ```
    /// ### Edge node
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(3);
    ///
    /// let node = momentum_lattice.get_node_by_index(&vec![9, 9, 9]);
    ///
    /// assert_eq!(node.get_neighbor_node(0).get_index(), &vec![9, 9, 9]);
    /// assert_eq!(node.get_neighbor_node(1).get_index(), &vec![0, 9, 9]);
    /// assert_eq!(node.get_neighbor_node(7).get_index(), &vec![0, 0, 9]);
    /// assert_eq!(node.get_neighbor_node(19).get_index(), &vec![0, 0, 0]);
    /// ```
    pub fn get_neighbor_node(&self, i: usize) -> Arc<Node> {
        self.get_neighbor_nodes()
            .get(&i)
            .cloned()
            .expect("Neighbor node not found")
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// let shallow_node = node.get_shallow_node();
    /// assert_eq!(shallow_node.get_density(), 1.0);
    /// assert_eq!(shallow_node.get_velocity(), vec![0.0, 0.0]);
    /// ```
    pub fn get_shallow_node(&self) -> &ShallowNode {
        &self.shallow_node
    }
}

impl Node {
    pub fn get_velocity_set_parameters(&self) -> &Arc<VelocitySetParameters> {
        &self.velocity_set_parameters
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_d(), &2);
    /// ```
    /// ## 3D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(3);
    ///
    /// assert_eq!(node.get_d(), &3);
    /// ```
    pub fn get_d(&self) -> &usize {
        &self.get_velocity_set_parameters().d
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_q(), &9);
    /// ```
    /// ## 3D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(3);
    ///
    /// assert_eq!(node.get_q(), &27);
    /// ```
    pub fn get_q(&self) -> &usize {
        &self.get_velocity_set_parameters().q
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// let c = node.get_c();
    ///
    /// assert_eq!(c[0], vec![0, 0]);
    /// assert_eq!(c[1], vec![1, 0]);
    /// assert_eq!(c[5], vec![1, 1]);
    /// ```
    /// ## 3D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(3);
    ///
    /// let c = node.get_c();
    ///
    /// assert_eq!(c[0], vec![0, 0, 0]);
    /// assert_eq!(c[1], vec![1, 0, 0]);
    /// assert_eq!(c[7], vec![1, 1, 0]);
    /// assert_eq!(c[19], vec![1, 1, 1]);
    /// ```
    pub fn get_c(&self) -> &Vec<Vec<i32>> {
        &self.get_velocity_set_parameters().c
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// let w = node.get_w();
    ///
    /// assert!(w[0] - 4.0 / 9.0 < 1e-12);
    /// assert!(w[1] - 1.0 / 9.0 < 1e-12);
    /// assert!(w[5] - 1.0 / 36.0 < 1e-12);
    /// ```
    /// ## 3D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(3);
    ///
    /// let w = node.get_w();
    ///
    /// assert!(w[0] - 8.0 / 27.0 < 1e-12);
    /// assert!(w[1] - 2.0 / 27.0 < 1e-12);
    /// assert!(w[7] - 1.0 / 54.0 < 1e-12);
    /// assert!(w[19] - 1.0 / 216.0 < 1e-12);
    /// ```
    pub fn get_w(&self) -> &Vec<Float> {
        &self.get_velocity_set_parameters().w
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// let q_bar = node.get_q_bar();
    ///
    /// assert_eq!(q_bar[0], 0);
    /// assert_eq!(q_bar[1], 3);
    /// assert_eq!(q_bar[5], 7);
    /// ```
    /// ## 3D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(3);
    ///
    /// let q_bar = node.get_q_bar();
    ///
    /// assert_eq!(q_bar[0], 0);
    /// assert_eq!(q_bar[1], 2);
    /// assert_eq!(q_bar[7], 8);
    /// assert_eq!(q_bar[19], 20);
    /// ```
    pub fn get_q_bar(&self) -> &Vec<usize> {
        &self.get_velocity_set_parameters().q_bar
    }

    pub fn get_q_faces(&self, boundary_face: &BoundaryFace) -> &Vec<usize> {
        self.get_velocity_set_parameters()
            .q_faces
            .get(boundary_face)
            .expect("Boundary face not found")
    }

    pub fn get_face_normal_direction(&self, boundary_face: &BoundaryFace) -> &usize {
        self.get_velocity_set_parameters()
            .face_normal_directions
            .get(boundary_face)
            .expect("Boundary face not found")
    }

    pub fn get_velocity_computation(&self) -> Option<VelocityComputation> {
        self.get_velocity_set_parameters().velocity_computation
    }

    pub fn get_opposite_direction(&self, direction: usize) -> usize {
        self.get_q_bar()[direction]
    }
}

impl Node {
    pub fn get_conversion_factor(&self) -> &Arc<ConversionFactor> {
        &self.conversion_factor
    }

    pub fn get_tau(&self) -> Float {
        self.get_conversion_factor().tau
    }

    pub fn get_physical_pressure(&self) -> Float {
        let density_prime = self.get_density() - LATTICE_DENSITY;
        let pressure_prime = CS_2 * density_prime;
        self.get_conversion_factor().reference_pressure
            + pressure_prime * self.get_conversion_factor().pressure_conversion_factor
    }

    pub fn get_physical_velocity(&self) -> Vec<Float> {
        self.get_velocity()
            .iter()
            .map(|u_x| u_x * self.get_conversion_factor().velocity_conversion_factor)
            .collect()
    }
}

impl Node {
    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_density(), 1.0);
    ///
    /// assert_eq!(node.get_f(), vec![0.0; 9]);
    /// node.compute_density();
    /// assert!((node.get_density() - 0.0).abs() < 1e-12);
    ///
    /// node.set_f(vec![0.1; 9]);
    /// assert_eq!(node.get_f(), vec![0.1; 9]);
    /// node.compute_density();
    /// assert!((node.get_density() - 0.9).abs() < 1e-12);
    ///
    /// node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    /// assert_eq!(node.get_f(), vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    /// node.compute_density();
    /// assert!((node.get_density() - 4.5).abs() < 1e-12);
    /// ```
    /// ## 3D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(3);
    ///
    /// assert_eq!(node.get_density(), 1.0);
    ///
    /// assert_eq!(node.get_f(), vec![0.0; 27]);
    /// node.compute_density();
    /// assert!((node.get_density() - 0.0).abs() < 1e-12);
    ///
    /// node.set_f(vec![0.1; 27]);
    /// assert_eq!(node.get_f(), vec![0.1; 27]);
    /// node.compute_density();
    /// assert!((node.get_density() - 2.7).abs() < 1e-12);
    /// ```
    pub fn compute_density(&self) {
        let f = self.get_f();
        let density = f.iter().sum::<Float>();
        self.set_density(density);
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_velocity(), vec![0.0; 2]);
    ///
    /// // Implicit computation:
    /// assert_eq!(node.get_f(), vec![0.0; 9]);
    /// node.compute_velocity(false);
    /// let actual = node.get_velocity();
    /// let target = vec![0.0; 2];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    /// node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    /// assert_eq!(node.get_f(), vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    /// node.compute_velocity(false);
    /// let actual = node.get_velocity();
    /// let target = vec![-0.2, -0.6];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    /// // Explicit computation:
    /// node.set_density(0.5);
    /// assert_eq!(node.get_density(), 0.5);
    /// node.compute_velocity(true);
    /// let actual = node.get_velocity();
    /// let target = vec![-0.4, -1.2];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    /// ```
    /// ## 3D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(3);
    ///
    /// assert_eq!(node.get_velocity(), vec![0.0; 3]);
    ///
    /// // Implicit computation:
    /// assert_eq!(node.get_f(), vec![0.0; 27]);
    /// node.compute_velocity(false);
    /// let actual = node.get_velocity();
    /// let target = vec![0.0; 3];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    /// node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]);
    /// assert_eq!(node.get_f(), vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]);
    /// node.compute_velocity(false);
    /// let actual = node.get_velocity();
    /// let target = vec![-0.1, 0.1, 0.3];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    ///
    /// // Explicit computation:
    /// node.set_density(0.5);
    /// assert_eq!(node.get_density(), 0.5);
    /// node.compute_velocity(true);
    /// let actual = node.get_velocity();
    /// let target = vec![-0.2, 0.2, 0.6];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    /// ```
    pub fn compute_velocity(&self, explicit_computation: bool) {
        let f = self.get_f();
        let density = self.get_density();
        let d = self.get_d();
        let c = self.get_c();
        match (explicit_computation, self.get_velocity_computation()) {
            (true, Some(velocity_computation)) => {
                let velocity = velocity_computation(density, f);
                self.set_velocity(velocity);
            }
            (_, _) => {
                let mut velocity = Vec::with_capacity(*d);
                (0..*d).for_each(|x| {
                    velocity.push(
                        f.iter()
                            .zip(c.iter())
                            .map(|(f_i, c_i)| f_i * ((*c_i)[x] as Float))
                            .sum::<Float>()
                            / density,
                    );
                });
                self.set_velocity(velocity);
            }
        }
    }

    /// $$ f\_{i}^{\text{eq}} = w\_{i}\rho\left[1+\frac{\mathbf{u}\cdot\mathbf{c}\_{i}}{c\_{s}^{2}}+\frac{\left(\mathbf{u}\cdot\mathbf{c}\_{i}\right)^{2}}{2 c\_{s}^{4}}-\frac{\mathbf{u}\cdot\mathbf{u}}{2 c\_{s}^{2}}\right] $$
    ///
    /// $$ f\_{0}^{\text{eq}} = \frac{2\rho}{9}\left(2-3\mathbf{u}^{2}\right) $$
    /// $$ f\_{1}^{\text{eq}} = \frac{\rho}{18}\left(2+6u\_{x}+9u\_{x}^{2}-3\mathbf{u}^{2}\right) $$
    /// $$ f\_{2}^{\text{eq}} = \frac{\rho}{18}\left(2+6u\_{y}+9u\_{y}^{2}-3\mathbf{u}^{2}\right) $$
    /// $$ f\_{3}^{\text{eq}} = \frac{\rho}{18}\left(2-6u\_{x}+9u\_{x}^{2}-3\mathbf{u}^{2}\right) $$
    /// $$ f\_{4}^{\text{eq}} = \frac{\rho}{18}\left(2-6u\_{y}+9u\_{y}^{2}-3\mathbf{u}^{2}\right) $$
    /// $$ f\_{5}^{\text{eq}} = \frac{\rho}{36}\left[1+3\left(u\_{x}+u\_{y}\right)+9u\_{x}u\_{y}+3\mathbf{u}^{2}\right] $$
    /// $$ f\_{6}^{\text{eq}} = \frac{\rho}{36}\left[1-3\left(u\_{x}-u\_{y}\right)-9u\_{x}u\_{y}+3\mathbf{u}^{2}\right] $$
    /// $$ f\_{7}^{\text{eq}} = \frac{\rho}{36}\left[1-3\left(u\_{x}+u\_{y}\right)+9u\_{x}u\_{y}+3\mathbf{u}^{2}\right] $$
    /// $$ f\_{8}^{\text{eq}} = \frac{\rho}{36}\left[1+3\left(u\_{x}-u\_{y}\right)-9u\_{x}u\_{y}+3\mathbf{u}^{2}\right] $$
    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// node.set_f(vec![0.1; 9]);
    /// node.compute_density(); // density = 0.9
    /// node.compute_equilibrium();
    /// let actual = node.get_f_eq();
    /// let target = vec![0.4, 0.1, 0.1, 0.1, 0.1, 0.025, 0.025, 0.025, 0.025];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    /// ```
    pub fn compute_equilibrium(&self) {
        let density = self.get_density();
        let velocity = self.get_velocity();
        let q = self.get_q();
        let c = self.get_c();
        let w = self.get_w();
        let mut f_eq = Vec::with_capacity(*q);
        let u_dot_u = velocity.iter().map(|u_x| u_x * u_x).sum::<Float>();
        (0..*q).for_each(|i| {
            let u_dot_c = velocity
                .iter()
                .zip(c[i].iter())
                .map(|(u_x, c_x)| u_x * (*c_x as Float))
                .sum::<Float>();
            f_eq.push(
                w[i] * density
                    * (1.0 + u_dot_c * CS_2_INV + 0.5 * u_dot_c * u_dot_c * CS_4_INV
                        - 0.5 * u_dot_u * CS_2_INV),
            );
        });
        self.set_f_eq(f_eq);
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::sync::Arc;
    /// use lbflow::constants::DELTA_T;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// node.set_f(vec![0.1; 9]);
    /// node.set_f_eq(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    ///
    /// node.compute_bgk_collision();
    /// let actual = node.get_f_star();
    /// let target = vec![0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    /// ```
    pub fn compute_bgk_collision(&self) {
        let tau = self.get_tau();
        let omega = DELTA_T / tau;
        let omega_prime = 1.0 - omega;
        let f = self.get_f();
        let f_eq = self.get_f_eq();
        let q = self.get_q();
        let mut f_star = Vec::with_capacity(*q);
        (0..*q).for_each(|i| {
            f_star.push(omega_prime * f[i] + omega * f_eq[i]);
        });
        self.set_f_star(f_star);
    }

    /// $$ f\_{i}(\mathbf{x}+\mathbf{c}\_{i}\Delta t, t+\Delta t) = f\_{i}^{\star}(\mathbf{x},t) $$
    /// # Examples
    /// ## Example 1
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(2);
    ///
    /// let node = momentum_lattice.get_node_by_index(&vec![3, 7]);
    ///
    /// node.get_neighbor_node(0).set_f_star(vec![0.1; 9]);
    /// node.get_neighbor_node(1).set_f_star(vec![0.2; 9]);
    /// node.get_neighbor_node(2).set_f_star(vec![0.3; 9]);
    /// node.get_neighbor_node(3).set_f_star(vec![0.4; 9]);
    /// node.get_neighbor_node(4).set_f_star(vec![0.5; 9]);
    /// node.get_neighbor_node(5).set_f_star(vec![0.6; 9]);
    /// node.get_neighbor_node(6).set_f_star(vec![0.7; 9]);
    /// node.get_neighbor_node(7).set_f_star(vec![0.8; 9]);
    /// node.get_neighbor_node(8).set_f_star(vec![0.9; 9]);
    ///
    /// node.compute_streaming();
    /// assert_eq!(node.get_f(), vec![0.1, 0.4, 0.5, 0.2, 0.3, 0.8, 0.9, 0.6, 0.7]);
    /// ```
    /// ## Example 2
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(2);
    ///
    /// let node = momentum_lattice.get_node_by_index(&vec![3, 7]);
    ///
    /// //                                          0    1    2    3    4    5    6    7    8
    /// node.get_neighbor_node(0).set_f_star(vec![1.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    /// node.get_neighbor_node(1).set_f_star(vec![0.1, 0.2, 0.3, 1.0, 0.5, 0.6, 0.7, 0.8, 0.9]);
    /// node.get_neighbor_node(2).set_f_star(vec![0.1, 0.2, 0.3, 0.4, 1.0, 0.6, 0.7, 0.8, 0.9]);
    /// node.get_neighbor_node(3).set_f_star(vec![0.1, 1.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    /// node.get_neighbor_node(4).set_f_star(vec![0.1, 0.2, 1.0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    /// node.get_neighbor_node(5).set_f_star(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 0.9]);
    /// node.get_neighbor_node(6).set_f_star(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]);
    /// node.get_neighbor_node(7).set_f_star(vec![0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 0.7, 0.8, 0.9]);
    /// node.get_neighbor_node(8).set_f_star(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0, 0.8, 0.9]);
    ///
    /// node.compute_streaming();
    /// assert_eq!(node.get_f(), vec![1.0; 9]);
    /// ```
    pub fn compute_streaming(&self) {
        let q = self.get_q();
        let mut f = vec![0.0; *q];
        self.get_neighbor_nodes()
            .iter()
            .for_each(|(i, neighbor_node)| {
                let i_bar = self.get_opposite_direction(*i);
                f[i_bar] = neighbor_node.get_f_star()[i_bar];
            });
        self.set_f(f);
    }

    pub fn compute_inner_bounce_back(&self) {
        let mut f = self.get_f();
        let f_star = self.get_f_star();
        self.get_bounce_back_neighbor_nodes()
            .iter()
            .for_each(|(&i, _)| {
                let i_bar = self.get_opposite_direction(i);
                f[i_bar] = f_star[i];
            });
        self.set_f(f);
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// node.set_density(0.95);
    /// node.set_velocity(vec![0.15, -0.25]);
    /// let (node_residual_density, node_residual_velocity) = node.compute_node_residuals();
    /// assert!((node_residual_density - (-0.05)).abs() < 1e-12);
    /// let actual = node_residual_velocity;
    /// let target = vec![0.15, -0.25];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    /// ```
    pub fn compute_node_residuals(&self) -> (Float, Vec<Float>) {
        let density = self.get_density();
        let velocity = self.get_velocity();
        let shallow_density = self.get_shallow_node().get_density();
        let shallow_velocity = self.get_shallow_node().get_velocity();
        let node_residual_density = density - shallow_density;
        let node_residual_velocity = velocity
            .iter()
            .zip(shallow_velocity.iter())
            .map(|(u_x, su_x)| (u_x - su_x))
            .collect::<Vec<Float>>();
        (node_residual_density, node_residual_velocity)
    }

    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// let shallow_node = node.get_shallow_node();
    /// assert_eq!(shallow_node.get_density(), 1.0);
    /// assert_eq!(shallow_node.get_velocity(), vec![0.0, 0.0]);
    ///
    /// node.set_density(0.9);
    /// node.set_velocity(vec![0.1, -0.1]);
    /// assert_eq!(shallow_node.get_density(), 1.0);
    /// assert_eq!(shallow_node.get_velocity(), vec![0.0, 0.0]);
    ///
    /// node.update_shallow_node();
    /// assert_eq!(shallow_node.get_density(), 0.9);
    /// assert_eq!(shallow_node.get_velocity(), vec![0.1, -0.1]);
    /// ```
    pub fn update_shallow_node(&self) {
        let density = self.get_density();
        let velocity = self.get_velocity();
        self.get_shallow_node().set_density(density);
        self.get_shallow_node().set_velocity(velocity);
    }
}

// ------------------------------------------------------------------- STRUCT: ShallowNode

#[derive(Debug)]
pub struct ShallowNode {
    pub density: RwLock<Float>,
    pub velocity: RwLock<Vec<Float>>,
}

impl ShallowNode {
    pub fn new(density: Float, velocity: Vec<Float>) -> Self {
        ShallowNode {
            density: RwLock::new(density),
            velocity: RwLock::new(velocity),
        }
    }
}

impl ShallowNode {
    /// # Examples
    /// ```
    /// # use lbflow::momentum::ShallowNode;
    /// let shallow_node = ShallowNode::new(1.0, vec![0.0, 0.0]);
    ///
    /// assert_eq!(shallow_node.get_density(), 1.0);
    /// ```
    pub fn get_density(&self) -> Float {
        // *self.density.borrow()
        *self.density.read().unwrap()
    }

    /// # Examples
    /// ```
    /// # use lbflow::momentum::ShallowNode;
    /// let shallow_node = ShallowNode::new(1.0, vec![0.0, 0.0]);
    ///
    /// shallow_node.set_density(0.9);
    /// assert_eq!(shallow_node.get_density(), 0.9);
    ///
    /// shallow_node.set_density(1.1);
    /// assert_eq!(shallow_node.get_density(), 1.1);
    /// ```
    pub fn set_density(&self, density: Float) {
        // self.density.replace(density);
        let mut density_guard = self.density.write().unwrap();
        *density_guard = density;
    }

    /// # Examples
    /// ```
    /// # use lbflow::momentum::ShallowNode;
    /// let shallow_node = ShallowNode::new(1.0, vec![0.0, 0.0]);
    ///
    /// assert_eq!(shallow_node.get_velocity(), vec![0.0, 0.0]);
    /// ```
    pub fn get_velocity(&self) -> Vec<Float> {
        // self.velocity.borrow().clone()
        self.velocity.read().unwrap().clone()
    }

    /// # Examples
    /// ```
    /// # use lbflow::momentum::ShallowNode;
    /// let shallow_node = ShallowNode::new(1.0, vec![0.0, 0.0]);
    ///
    /// shallow_node.set_velocity(vec![0.1, -0.1]);
    /// assert_eq!(shallow_node.get_velocity(), vec![0.1, -0.1]);
    ///
    /// shallow_node.set_velocity(vec![0.3, 0.7]);
    /// assert_eq!(shallow_node.get_velocity(), vec![0.3, 0.7]);
    /// ```
    pub fn set_velocity(&self, velocity: Vec<Float>) {
        // self.velocity.replace(velocity);
        let mut velocity_guard = self.velocity.write().unwrap();
        *velocity_guard = velocity;
    }
}
