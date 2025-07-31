use super::constants::*;
use super::momentum::Lattice as MomentumLattice;
use super::momentum::Node as MomentumNode;
use super::{BoundaryFace, NodeType};
use std::cell::RefCell;
use std::rc::Rc;

#[derive(Debug)]
pub struct Node {
    pub momentum_node: Rc<MomentumNode>,
    pub temperature: RefCell<Float>,
    pub g: RefCell<Vec<Float>>,
    pub g_eq: RefCell<Vec<Float>>,
    pub g_star: RefCell<Vec<Float>>,
}

#[derive(Debug)]
pub struct Lattice {
    pub momentum_lattice: Rc<MomentumLattice>,
    pub nodes: Vec<Rc<Node>>,
    pub tau_g: Float,
}

impl Node {
    pub fn new(momentum_node: Rc<MomentumNode>, temperature: Float) -> Self {
        let q_par = *momentum_node.get_q();
        Node {
            momentum_node,
            temperature: RefCell::new(temperature),
            g: RefCell::new(vec![0.0; q_par]),
            g_eq: RefCell::new(vec![0.0; q_par]),
            g_star: RefCell::new(vec![0.0; q_par]),
        }
    }
}

impl Node {
    pub fn get_momentum_node(&self) -> &Rc<MomentumNode> {
        &self.momentum_node
    }

    pub fn get_density(&self) -> Float {
        self.get_momentum_node().get_density()
    }

    pub fn get_velocity(&self) -> Vec<Float> {
        self.get_momentum_node().get_velocity()
    }

    pub fn get_temperature(&self) -> Float {
        self.temperature.borrow().clone()
    }

    pub fn set_temperature(&self, temperature: Float) {
        self.temperature.replace(temperature);
    }

    pub fn get_g(&self) -> Vec<Float> {
        self.g.borrow().clone()
    }

    fn set_g(&self, g: Vec<Float>) {
        self.g.replace(g);
    }

    pub fn get_g_eq(&self) -> Vec<Float> {
        self.g_eq.borrow().clone()
    }

    fn set_g_eq(&self, g_eq: Vec<Float>) {
        self.g_eq.replace(g_eq);
    }

    pub fn get_g_star(&self) -> Vec<Float> {
        self.g_star.borrow().clone()
    }

    fn set_g_star(&self, g_star: Vec<Float>) {
        self.g_star.replace(g_star);
    }
}

impl Node {
    pub fn get_d(&self) -> &usize {
        &self.get_momentum_node().get_velocity_set_parameters().d
    }

    pub fn get_q(&self) -> &usize {
        &self.get_momentum_node().get_velocity_set_parameters().q
    }

    pub fn get_c(&self) -> &Vec<Vec<i32>> {
        &self.get_momentum_node().get_velocity_set_parameters().c
    }

    pub fn get_w(&self) -> &Vec<Float> {
        &self.get_momentum_node().get_velocity_set_parameters().w
    }

    pub fn get_q_bar(&self) -> &Vec<usize> {
        &self.get_momentum_node().get_velocity_set_parameters().q_bar
    }
}

impl Lattice {
    pub fn new(momentum_lattice: Rc<MomentumLattice>, tau_g: Float) -> Self {
        let nodes = momentum_lattice
            .get_nodes()
            .iter()
            .map(|node| {
                Rc::new(Node::new(Rc::clone(node), 300.0)) // Default temperature
            })
            .collect();

        Lattice {
            momentum_lattice,
            nodes,
            tau_g,
        }
    }
}
