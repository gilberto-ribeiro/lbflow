// ------------------------------------------------------------------------------- MODULES

pub mod bc;
pub mod io;
mod legacy;

// ------------------------------------------------------------------------------- IMPORTS

use crate::parameters::{self, *};
use crate::simulation::Simulation;

use super::constants::*;
use super::{BoundaryFace, NodeType};
use super::{FACES_2D, FACES_3D};
use super::{VelocitySet, VelocitySetParameters};
use bc::BoundaryCondition;
use core::num;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

// -------------------------------------------------------------------- STRUCT: Parameters

pub struct Parameters {
    pub n: Vec<usize>,
    pub tau: Float,
    pub physical_delta_x: Float,
    pub physical_delta_t: Float,
    pub velocity_set: VelocitySet,
    pub boundary_conditions: Vec<(BoundaryFace, BoundaryCondition)>,
    pub node_types: Vec<NodeType>,
}

impl Default for Parameters {
    fn default() -> Self {
        Parameters {
            n: vec![10, 10],
            tau: 0.5,
            physical_delta_x: 0.01,
            physical_delta_t: 0.01,
            velocity_set: VelocitySet::D2Q9,
            node_types: parameters::functions::only_fluid_nodes(vec![10, 10]),
            boundary_conditions: vec![
                (BoundaryFace::West, BoundaryCondition::NoSlip),
                (BoundaryFace::East, BoundaryCondition::NoSlip),
                (BoundaryFace::South, BoundaryCondition::NoSlip),
                (
                    BoundaryFace::North,
                    BoundaryCondition::BounceBack {
                        density: 1.0,
                        velocity: vec![0.1, 0.0],
                    },
                ),
            ],
        }
    }
}

impl Parameters {
    pub fn test_default(dim: usize) -> Self {
        match dim {
            2 => Default::default(),
            3 => Parameters {
                n: vec![10, 10, 10],
                velocity_set: VelocitySet::D3Q27,
                node_types: parameters::functions::only_fluid_nodes(vec![10, 10, 10]),
                boundary_conditions: vec![
                    (BoundaryFace::West, BoundaryCondition::NoSlip),
                    (BoundaryFace::East, BoundaryCondition::NoSlip),
                    (BoundaryFace::South, BoundaryCondition::NoSlip),
                    (
                        BoundaryFace::North,
                        BoundaryCondition::BounceBack {
                            density: 1.0,
                            velocity: vec![0.1, 0.0, 0.0],
                        },
                    ),
                    (BoundaryFace::Bottom, BoundaryCondition::NoSlip),
                    (BoundaryFace::Top, BoundaryCondition::NoSlip),
                ],
                ..Default::default()
            },
            _ => panic!("Unsupported dimension: {dim}"),
        }
    }
}

// --------------------------------------------------------------------- STRUCT: Residuals

#[derive(Debug, Clone)]
pub struct Residuals {
    pub density: Float,
    pub velocity: Vec<Float>,
}

impl Residuals {
    pub fn new(density: Float, velocity: Vec<Float>) -> Self {
        Residuals { density, velocity }
    }
}

impl Residuals {
    pub fn get_density(&self) -> Float {
        self.density
    }

    pub fn set_density(&mut self, density: Float) {
        self.density = density;
    }

    pub fn get_velocity(&self) -> &Vec<Float> {
        &self.velocity
    }

    pub fn set_velocity(&mut self, velocity: Vec<Float>) {
        self.velocity = velocity;
    }
}

// -------------------------------------------------------------------------- STRUCT: Node

#[derive(Debug, PartialEq)]
pub struct Node {
    density: RefCell<Float>,
    velocity: RefCell<Vec<Float>>,
    f: RefCell<Vec<Float>>,
    f_eq: RefCell<Vec<Float>>,
    f_star: RefCell<Vec<Float>>,
    node_type: NodeType,
    index: Vec<usize>,
    coordinates: Vec<Float>,
    velocity_set_parameters: Rc<VelocitySetParameters>,
    neighbor_nodes: RefCell<Option<HashMap<usize, Rc<Node>>>>,
    bounce_back_neighbor_nodes: RefCell<Option<HashMap<usize, Rc<Node>>>>,
    shallow_node: ShallowNode,
}

impl Node {
    pub fn new(
        density: Float,
        velocity: Vec<Float>,
        node_type: NodeType,
        index: Vec<usize>,
        coordinates: Vec<Float>,
        velocity_set_parameters: Rc<VelocitySetParameters>,
    ) -> Self {
        let q = velocity_set_parameters.q;
        Node {
            density: RefCell::new(density),
            velocity: RefCell::new(velocity.clone()),
            f: RefCell::new(vec![0.0; q]),
            f_eq: RefCell::new(vec![0.0; q]),
            f_star: RefCell::new(vec![0.0; q]),
            node_type,
            index,
            coordinates,
            velocity_set_parameters,
            neighbor_nodes: RefCell::new(None),
            bounce_back_neighbor_nodes: RefCell::new(None),
            shallow_node: ShallowNode::new(density, velocity.clone()),
        }
    }

    pub fn test_default(dim: usize) -> Self {
        match dim {
            2 => Node::new(
                1.0,
                vec![0.0, 0.0],
                NodeType::Fluid,
                vec![3, 7],
                vec![0.035, 0.075],
                Rc::new(VelocitySetParameters::test_default(2)),
            ),
            3 => Node::new(
                1.0,
                vec![0.0, 0.0, 0.0],
                NodeType::Fluid,
                vec![3, 5, 7],
                vec![0.035, 0.055, 0.075],
                Rc::new(VelocitySetParameters::test_default(3)),
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
            Rc::new(VelocitySetParameters::default()),
        )
    }
}

impl Node {
    /// # Examples
    /// ```
    /// # use std::rc::Rc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_density(), 1.0);
    /// assert_eq!(&node.get_density(), &1.0);
    /// ```
    pub fn get_density(&self) -> Float {
        *self.density.borrow()
    }

    /// # Examples
    /// ```
    /// # use std::rc::Rc;
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
    pub fn set_density(&self, density: Float) {
        self.density.replace(density);
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
        self.velocity.borrow().clone()
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    pub fn set_velocity(&self, velocity: Vec<Float>) {
        self.velocity.replace(velocity);
    }

    /// # Examples
    /// ```
    /// # use std::rc::Rc;
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
        self.f.borrow().clone()
    }

    /// # Examples
    /// ```
    /// # use std::rc::Rc;
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
        self.f.replace(f);
    }

    /// # Examples
    /// ```
    /// # use std::rc::Rc;
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
        self.f_eq.borrow().clone()
    }

    /// # Examples
    /// ```
    /// # use std::rc::Rc;
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
    pub fn set_f_eq(&self, f_eq: Vec<Float>) {
        self.f_eq.replace(f_eq);
    }

    /// # Examples
    /// ```
    /// # use std::rc::Rc;
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
        self.f_star.borrow().clone()
    }

    /// # Examples
    /// ```
    /// # use std::rc::Rc;
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
    pub fn set_f_star(&self, f_star: Vec<Float>) {
        self.f_star.replace(f_star);
    }

    /// # Examples
    /// ```
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    pub fn get_neighbor_nodes(&self) -> HashMap<usize, Rc<Node>> {
        self.neighbor_nodes.borrow().as_ref().cloned().unwrap()
    }

    pub fn get_bounce_back_neighbor_nodes(&self) -> HashMap<usize, Rc<Node>> {
        self.bounce_back_neighbor_nodes
            .borrow()
            .as_ref()
            .cloned()
            .unwrap()
    }

    /// # Examples
    /// ## 2D
    /// ### Inner node
    /// ```
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    pub fn get_neighbor_node(&self, i: usize) -> Rc<Node> {
        self.get_neighbor_nodes()
            .get(&i)
            .cloned()
            .expect("Neighbor node not found")
    }

    /// # Examples
    /// ```
    /// # use std::rc::Rc;
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
    pub fn get_velocity_set_parameters(&self) -> &Rc<VelocitySetParameters> {
        &self.velocity_set_parameters
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::rc::Rc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_d(), &2);
    /// ```
    /// ## 3D
    /// ```
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_q(), &9);
    /// ```
    /// ## 3D
    /// ```
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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

    pub fn get_velocity_computation(&self) -> Option<fn(Float, Vec<Float>) -> Vec<Float>> {
        self.get_velocity_set_parameters().velocity_computation
    }

    pub fn get_opposite_direction(&self, direction: usize) -> usize {
        self.get_q_bar()[direction]
    }
}

impl Node {
    /// # Examples
    /// ## 2D
    /// ```
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
        (0..*q).for_each(|i| {
            let u_dot_c = velocity
                .iter()
                .zip(c[i].iter())
                .map(|(u_x, c_x)| u_x * (*c_x as Float))
                .sum::<Float>();
            let u_dot_u = velocity.iter().map(|u_x| u_x * u_x).sum::<Float>();
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
    /// # use std::rc::Rc;
    /// use lbflow::constants::DELTA_T;
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// # use lbflow::velocity_set::VelocitySetParameters;
    /// let node = Node::test_default(2);
    ///
    /// node.set_f(vec![0.1; 9]);
    /// node.set_f_eq(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
    ///
    /// let tau = 0.5;
    /// let omega = DELTA_T / tau;
    /// let omega_prime = 1.0 - omega;
    ///
    /// node.compute_bgk_collision(omega, omega_prime);
    /// let actual = node.get_f_star();
    /// let target = vec![0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7];
    /// for (a, b) in actual.iter().zip(target.iter()) {
    ///     assert!((a - b).abs() < 1e-12);
    /// }
    /// ```
    pub fn compute_bgk_collision(&self, omega: Float, omega_prime: Float) {
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
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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
            .for_each(|(i, _)| {
                let i_bar = self.get_opposite_direction(*i);
                f[i_bar] = f_star[*i];
            });
        self.set_f(f);
    }

    /// # Examples
    /// ```
    /// # use std::rc::Rc;
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
    /// # use std::rc::Rc;
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

#[derive(Debug, PartialEq, Clone)]
pub struct ShallowNode {
    pub density: RefCell<Float>,
    pub velocity: RefCell<Vec<Float>>,
}

impl ShallowNode {
    pub fn new(density: Float, velocity: Vec<Float>) -> Self {
        ShallowNode {
            density: RefCell::new(density),
            velocity: RefCell::new(velocity),
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
        *self.density.borrow()
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
        self.density.replace(density);
    }

    /// # Examples
    /// ```
    /// # use lbflow::momentum::ShallowNode;
    /// let shallow_node = ShallowNode::new(1.0, vec![0.0, 0.0]);
    ///
    /// assert_eq!(shallow_node.get_velocity(), vec![0.0, 0.0]);
    /// ```
    pub fn get_velocity(&self) -> Vec<Float> {
        self.velocity.borrow().clone()
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
        self.velocity.replace(velocity);
    }
}

// ----------------------------------------------------------------------- STRUCT: Lattice

#[derive(Debug)]
pub struct Lattice {
    nodes: Vec<Rc<Node>>,
    pub n: Vec<usize>,
    pub physical_delta_x: Float,
    pub physical_delta_t: Float,
    pub tau: Float,
    pub velocity_set_parameters: Rc<VelocitySetParameters>,
    fluid_nodes: Vec<Rc<Node>>,
    boundary_nodes: HashMap<BoundaryFace, Vec<Rc<Node>>>,
    bounce_back_nodes: Vec<Rc<Node>>,
    pub boundary_conditions: HashMap<BoundaryFace, BoundaryCondition>,
    residuals: RefCell<Residuals>,
    time_step: RefCell<usize>,
}

impl Lattice {
    pub fn new(parameters: Parameters) -> Self {
        println!(
            "Selecting velocity set for the lattice: {:?}\n",
            parameters.velocity_set
        );
        let velocity_set_parameters =
            Rc::new(super::get_velocity_set_parameters(parameters.velocity_set));

        let d = &velocity_set_parameters.d;
        let c = &(velocity_set_parameters.c);
        let n = parameters.n;
        let node_types = parameters.node_types;
        let physical_delta_x = parameters.physical_delta_x;
        let num_nodes = n.iter().product::<usize>();
        if num_nodes != node_types.len() {
            panic!(
                "Number of nodes ({num_nodes}) does not match the length of node types ({})",
                node_types.len()
            );
        }
        let (boundary_faces, limit_index, dimension_index) = match d {
            2 => (
                FACES_2D.to_vec(),
                vec![0, n[0] - 1, 0, n[1] - 1],
                vec![0, 0, 1, 1],
            ),
            3 => (
                FACES_3D.to_vec(),
                vec![0, n[0] - 1, 0, n[1] - 1, 0, n[2] - 1],
                vec![0, 0, 1, 1, 2, 2],
            ),
            _ => panic!("Unsupported dimension: {d}"),
        };

        println!("Creating lattice with dimensions: {n:?}\n");
        let nodes = (0..num_nodes)
            .zip(node_types.iter())
            .map(|(i, &node_type)| {
                let index = match d {
                    2 => {
                        let x = i % n[0];
                        let y = i / n[0];
                        vec![x, y]
                    }
                    3 => {
                        let x = i % n[0];
                        let y = (i / n[0]) % n[1];
                        let z = i / (n[0] * n[1]);
                        vec![x, y, z]
                    }
                    _ => panic!("Unsupported dimension: {d}"),
                };
                let coordinates = index
                    .iter()
                    .map(|&x| (x as Float + 0.5) * physical_delta_x)
                    .collect::<Vec<Float>>();
                Rc::new(Node::new(
                    1.0,
                    vec![0.0; *d],
                    node_type,
                    index,
                    coordinates,
                    Rc::clone(&velocity_set_parameters),
                ))
            })
            .collect::<Vec<Rc<Node>>>();

        let nx = n[0];
        let ny = n[1];
        let nz = match d {
            2 => 1,
            3 => n[2],
            _ => panic!("Unsupported dimension: {d}"),
        };

        let mut nodes_matrix: Vec<Vec<Vec<Option<Rc<Node>>>>> = vec![vec![vec![None; nx]; ny]; nz];
        nodes.iter().for_each(|node| {
            let index = node.get_index();
            let x = index[0];
            let y = index[1];
            let z = match d {
                2 => 0,
                3 => index[2],
                _ => panic!("Unsupported dimension: {d}"),
            };
            nodes_matrix[z][y][x] = Some(Rc::clone(node));
        });

        println!("Setting up neighbor nodes for each node in the lattice...\n");
        nodes.iter().enumerate().for_each(|(n_i, node)| {
            crate::io::progress_bar(n_i, num_nodes);
            let index = node.get_index();
            let neighbors = c
                .iter()
                .enumerate()
                .map(|(i, c_i)| {
                    let neighbor_index = index
                        .iter()
                        .zip(c_i.iter().zip(n.iter()))
                        .map(|(idx, (&c_x, &n_x))| {
                            (*idx as i32 + c_x).rem_euclid(n_x as i32) as usize
                        })
                        .collect::<Vec<usize>>();
                    // let neighbor_node = nodes
                    //     .iter()
                    //     .find(|node| node.get_index() == &neighbor_index)
                    //     .cloned()
                    //     .unwrap();
                    let x = neighbor_index[0];
                    let y = neighbor_index[1];
                    let z = match d {
                        2 => 0,
                        3 => neighbor_index[2],
                        _ => panic!("Unsupported dimension: {d}"),
                    };
                    let neighbor_node = Rc::clone(nodes_matrix[z][y][x].as_ref().unwrap());
                    (i, neighbor_node)
                })
                .collect::<HashMap<usize, Rc<Node>>>();
            node.neighbor_nodes.replace(Some(neighbors));
        });

        println!("Identifying fluid nodes...\n");
        let fluid_nodes = nodes
            .iter()
            .filter(|node| matches!(node.get_node_type(), NodeType::Fluid))
            .cloned()
            .collect::<Vec<Rc<Node>>>();

        println!("Identifying boundary nodes...\n");
        let boundary_nodes = boundary_faces
            .iter()
            .zip(limit_index.iter().zip(dimension_index.iter()))
            .map(|(face, (&limit, &x))| {
                let face_nodes = fluid_nodes
                    .iter()
                    .filter(|node| node.get_index()[x] == limit)
                    .cloned()
                    .collect::<Vec<Rc<Node>>>();
                (*face, face_nodes)
            })
            .collect::<HashMap<BoundaryFace, Vec<Rc<Node>>>>();

        println!("Identifying bounce-back nodes...\n");
        let mut bounce_back_nodes = Vec::new();
        fluid_nodes.iter().for_each(|node| {
            let bounce_back_neighbor_nodes = node
                .get_neighbor_nodes()
                .iter()
                .filter(|(_, neighbor_node)| {
                    matches!(neighbor_node.get_node_type(), NodeType::Solid)
                })
                .map(|(i, neighbor_node)| (*i, Rc::clone(neighbor_node)))
                .collect::<HashMap<usize, Rc<Node>>>();
            if !bounce_back_neighbor_nodes.is_empty() {
                node.bounce_back_neighbor_nodes
                    .replace(Some(bounce_back_neighbor_nodes));
                bounce_back_nodes.push(Rc::clone(node));
            }
        });

        Lattice {
            nodes,
            n,
            physical_delta_x: parameters.physical_delta_x,
            physical_delta_t: parameters.physical_delta_t,
            tau: parameters.tau,
            velocity_set_parameters: Rc::clone(&velocity_set_parameters),
            fluid_nodes,
            boundary_nodes,
            boundary_conditions: HashMap::from_iter(parameters.boundary_conditions),
            bounce_back_nodes,
            residuals: RefCell::new(Residuals::new(0.0, vec![0.0; *d])),
            time_step: RefCell::new(0),
        }
    }

    pub fn test_default(dim: usize) -> Self {
        let parameters = Parameters::test_default(dim);
        Lattice::new(parameters)
    }
}

impl Default for Lattice {
    fn default() -> Self {
        let parameters = Parameters::default();
        Lattice::new(parameters)
    }
}

impl Lattice {
    /// # Examples
    /// ```
    /// # use std::rc::Rc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(2);
    ///
    /// assert_eq!(momentum_lattice.get_node_by_index(&vec![3, 7]).get_index(), &vec![3, 7]);
    pub fn get_node_by_index(&self, index: &Vec<usize>) -> &Rc<Node> {
        self.nodes
            .iter()
            .find(|node| node.get_index() == index)
            .unwrap()
    }

    /// # Examples
    /// ```
    /// # use lbflow::momentum::Lattice;
    /// # use lbflow::BoundaryFace;
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(2);
    ///
    /// let west_nodes = momentum_lattice.get_boundary_nodes_by_face(&BoundaryFace::West)
    ///     .iter()
    ///     .map(|node| node.get_index())
    ///     .collect::<Vec<&Vec<usize>>>();
    /// let east_nodes = momentum_lattice.get_boundary_nodes_by_face(&BoundaryFace::East)
    ///     .iter()
    ///     .map(|node| node.get_index())
    ///     .collect::<Vec<&Vec<usize>>>();
    /// let south_nodes = momentum_lattice.get_boundary_nodes_by_face(&BoundaryFace::South)
    ///     .iter()
    ///     .map(|node| node.get_index())
    ///     .collect::<Vec<&Vec<usize>>>();
    /// let north_nodes = momentum_lattice.get_boundary_nodes_by_face(&BoundaryFace::North)
    ///     .iter()
    ///     .map(|node| node.get_index())
    ///     .collect::<Vec<&Vec<usize>>>();
    ///
    /// assert!(west_nodes.contains(&&vec![0, 0]));
    /// assert!(west_nodes.contains(&&vec![0, 5]));
    /// assert!(west_nodes.contains(&&vec![0, 9]));
    ///
    /// assert!(east_nodes.contains(&&vec![9, 0]));
    /// assert!(east_nodes.contains(&&vec![9, 5]));
    /// assert!(east_nodes.contains(&&vec![9, 9]));
    ///
    /// assert!(south_nodes.contains(&&vec![0, 0]));
    /// assert!(south_nodes.contains(&&vec![5, 0]));
    /// assert!(south_nodes.contains(&&vec![9, 0]));
    ///
    /// assert!(north_nodes.contains(&&vec![0, 9]));
    /// assert!(north_nodes.contains(&&vec![5, 9]));
    /// assert!(north_nodes.contains(&&vec![9, 9]));
    ///
    /// let node = momentum_lattice.get_node_by_index(&vec![3, 7]);
    ///
    /// assert!(!west_nodes.contains(&node.get_index()));
    /// assert!(!east_nodes.contains(&node.get_index()));
    /// assert!(!south_nodes.contains(&node.get_index()));
    /// assert!(!north_nodes.contains(&node.get_index()));
    ///
    /// let neighbor_node = node.get_neighbor_node(2).get_neighbor_node(2);
    /// assert_eq!(neighbor_node.get_index(), &vec![3, 9]);
    /// assert!(north_nodes.contains(&neighbor_node.get_index()));
    /// ```
    pub fn get_boundary_nodes_by_face(&self, boundary_face: &BoundaryFace) -> &Vec<Rc<Node>> {
        self.boundary_nodes
            .get(boundary_face)
            .expect("Boundary face not found")
    }

    /// # Examples
    /// ```
    /// # use lbflow::momentum::Lattice;
    /// # use lbflow::BoundaryFace;
    /// # use lbflow::momentum::bc::BoundaryCondition;
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(2);
    ///
    /// assert_eq!(momentum_lattice.get_boundary_condition(&BoundaryFace::West), &BoundaryCondition::NoSlip);
    /// assert_eq!(momentum_lattice.get_boundary_condition(&BoundaryFace::East), &BoundaryCondition::NoSlip);
    /// assert_eq!(momentum_lattice.get_boundary_condition(&BoundaryFace::South), &BoundaryCondition::NoSlip);
    /// assert_eq!(momentum_lattice.get_boundary_condition(&BoundaryFace::North),
    ///     &BoundaryCondition::BounceBack {
    ///         density: 1.0,
    ///         velocity: vec![0.1, 0.0],
    ///     }
    /// );
    /// ```
    pub fn get_boundary_condition(&self, boundary_face: &BoundaryFace) -> &BoundaryCondition {
        self.boundary_conditions
            .get(boundary_face)
            .expect("Boundary faces not found")
    }

    pub fn get_nodes(&self) -> &Vec<Rc<Node>> {
        &self.nodes
    }

    /// Examples
    /// 2D
    /// ```
    /// # use std::rc::Rc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(2);
    ///
    /// assert_eq!(momentum_lattice.get_node(0).get_index(), &vec![0, 0]);
    /// assert_eq!(momentum_lattice.get_node(1).get_index(), &vec![1, 0]);
    /// assert_eq!(momentum_lattice.get_node(73).get_index(), &vec![3, 7]);
    ///
    /// let node_0 = momentum_lattice.get_node(73);
    /// let node_1 = momentum_lattice.get_node(73);
    ///
    /// assert_eq!(node_0.get_index(), &vec![3, 7]);
    /// assert_eq!(node_1.get_index(), &vec![3, 7]);
    /// ```
    pub fn get_node(&self, i: usize) -> &Rc<Node> {
        &self.nodes[i]
    }

    pub fn get_fluid_nodes(&self) -> &Vec<Rc<Node>> {
        &self.fluid_nodes
    }

    pub fn get_bounce_back_nodes(&self) -> &Vec<Rc<Node>> {
        &self.bounce_back_nodes
    }

    pub fn get_boundary_nodes(&self) -> &HashMap<BoundaryFace, Vec<Rc<Node>>> {
        &self.boundary_nodes
    }

    pub fn get_time_step(&self) -> usize {
        *self.time_step.borrow()
    }

    pub fn get_residuals(&self) -> Residuals {
        self.residuals.borrow().clone()
    }

    pub fn set_residuals(&self, residuals: Residuals) {
        self.residuals.replace(residuals);
    }

    pub fn get_nx(&self) -> usize {
        self.n[0]
    }

    pub fn get_ny(&self) -> usize {
        self.n[1]
    }

    pub fn get_nz(&self) -> usize {
        match self.get_d() {
            2 => self.n[0],
            3 => self.n[1],
            _ => panic!("Unsupported dimension: {}", self.get_d()),
        }
    }

    pub fn get_number_of_nodes(&self) -> usize {
        self.n.iter().product()
    }
}

impl Lattice {
    pub fn get_d(&self) -> &usize {
        &self.velocity_set_parameters.d
    }
}

impl Lattice {
    pub fn initialize_nodes(&self) {
        self.get_fluid_nodes().iter().for_each(|node| {
            node.compute_equilibrium();
            node.set_f(node.get_f_eq());
        });
    }

    pub fn update_density_and_velocity_step(&self) {
        self.get_fluid_nodes().iter().for_each(|node| {
            node.compute_density();
            node.compute_velocity(false);
        });
    }

    pub fn equilibrium_step(&self) {
        self.get_fluid_nodes().iter().for_each(|node| {
            node.compute_equilibrium();
        });
    }

    pub fn bgk_collision_step(&self) {
        let omega = DELTA_T / self.tau;
        let omega_prime = 1.0 - omega;
        self.get_fluid_nodes().iter().for_each(|node| {
            node.compute_bgk_collision(omega, omega_prime);
        });
    }

    pub fn streaming_step(&self) {
        self.get_fluid_nodes().iter().for_each(|node| {
            node.compute_streaming();
        });
    }

    pub fn inner_bounce_back_step(&self) {
        self.get_bounce_back_nodes().iter().for_each(|node| {
            node.compute_inner_bounce_back();
        });
    }

    pub fn boundary_conditions_step(&self) {
        self.get_boundary_nodes()
            .iter()
            .for_each(|(boundary_face, nodes)| {
                let boundary_condition = self.get_boundary_condition(boundary_face);
                match boundary_condition {
                    BoundaryCondition::NoSlip => {
                        nodes.iter().for_each(|node| {
                            node.compute_no_slip_bc(boundary_face);
                        });
                    }
                    BoundaryCondition::BounceBack { density, velocity } => {
                        nodes.iter().for_each(|node| {
                            node.compute_bounce_back_bc(boundary_face, density, velocity);
                        });
                    }
                    BoundaryCondition::AntiBounceBack { density } => {
                        nodes.iter().for_each(|node| {
                            node.compute_anti_bounce_back_bc(boundary_face, density);
                        });
                    }
                    BoundaryCondition::Periodic => {}
                }
            });
    }

    pub fn compute_lattice_residuals(&self) {
        let d = self.get_d();
        let node_residuals = self
            .get_fluid_nodes()
            .iter()
            .map(|node| node.compute_node_residuals())
            .collect::<Vec<(Float, Vec<Float>)>>();
        let residuals_density = node_residuals
            .iter()
            .map(|(residual_density, _)| residual_density * residual_density)
            .sum::<Float>()
            .sqrt();
        let residuals_velocity = (0..*d)
            .map(|x| {
                node_residuals
                    .iter()
                    .map(|(_, residual_velocity)| residual_velocity[x] * residual_velocity[x])
                    .sum::<Float>()
                    .sqrt()
            })
            .collect::<Vec<Float>>();
        self.set_residuals(Residuals {
            density: residuals_density,
            velocity: residuals_velocity,
        });
    }

    pub fn update_shallow_nodes(&self) {
        self.get_fluid_nodes().iter().for_each(|node| {
            node.update_shallow_node();
        });
    }

    pub fn next_time_step(&self) {
        *self.time_step.borrow_mut() += 1;
    }
}

// ----------------------------------------------------------------------------- FUNCTIONS

pub fn run() {
    let _ = crate::io::create_case_directories();
    let simulation = Simulation::new(get_simulation_parameters());
    let lattice = Lattice::new(get_momentum_parameters());
    let _ = lattice.write_coordinates();
    lattice.initialize_nodes();
    loop {
        lattice.update_density_and_velocity_step();
        lattice.equilibrium_step();
        lattice.bgk_collision_step();
        lattice.streaming_step();
        lattice.inner_bounce_back_step();
        lattice.boundary_conditions_step();
        lattice.write_data(simulation.get_data_write_mode());
        lattice.compute_lattice_residuals();
        lattice.print_residuals();
        let _ = lattice.write_residuals();
        lattice.update_shallow_nodes();
        lattice.next_time_step();
    }
}
