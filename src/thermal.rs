use std::rc::Rc;

use super::Float;
use super::momentum::Node as MomentumNode;
use super::momentum::Lattice as MomentumLattice;

#[derive(Debug)]
pub struct Node {
    pub momentum_node: Rc<MomentumNode>,
    pub temperature: Float,
    pub g: [Float; 9],
    pub g_eq: [Float; 9],
    pub g_star: [Float; 9],
}

#[derive(Debug)]
pub struct Lattice {
    pub momentum_lattice: Rc<MomentumLattice>,
    pub nodes: Vec<Rc<Node>>,
    pub tau_g: Float,
}

impl Node {
    pub fn new(momentum_node: Rc<MomentumNode>, temperature: Float) -> Self {
        Node {
            momentum_node,
            temperature,
            g: [0.0; 9],
            g_eq: [0.0; 9],
            g_star: [0.0; 9],
        }
    }
}

impl Lattice {
    pub fn new(momentum_lattice: Rc<MomentumLattice>, tau_g: Float) -> Self {
        let nodes = momentum_lattice.nodes.iter().map(|node| {
            Rc::new(Node::new(Rc::clone(node), 300.0)) // Default temperature
        }).collect();
        
        Lattice {
            momentum_lattice,
            nodes,
            tau_g,
        }
    }
}