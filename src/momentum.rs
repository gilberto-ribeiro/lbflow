use std::rc::Rc;

use super::Float;

#[derive(Debug)]
pub struct Node {
    pub density: Float,
    pub velocity: [Float; 2],
    pub f: [Float; 9],
    pub f_eq: [Float; 9],
    pub f_star: [Float; 9],
    pub node_type: NodeType,
    pub index: [usize; 2],
}

#[derive(Debug)]
pub struct Lattice {
    pub nodes: Vec<Rc<Node>>,
    pub nx: usize,
    pub ny: usize,
    pub delta_x: Float,
    pub delta_t: Float,
    pub tau: Float,
}

#[derive(Debug)]
pub enum NodeType {
    Fluid,
    Solid,
}

impl Node {
    pub fn new(density: Float, velocity: [Float; 2], node_type: NodeType, index: [usize; 2]) -> Self {
        Node {
            density,
            velocity,
            f: [0.0; 9],
            f_eq: [0.0; 9],
            f_star: [0.0; 9],
            node_type,
            index,
        }
    }
}

impl Lattice {
    pub fn new(nx: usize, ny: usize, delta_x: Float, delta_t: Float, tau: Float) -> Self {
        let nodes = (0..nx * ny)
            .map(|i| {
                let x = i % nx;
                let y = i / nx;
                Rc::new(Node::new(1.0, [0.0, 0.0], NodeType::Fluid, [x, y]))
            })
            .collect();
        
        Lattice {
            nodes,
            nx,
            ny,
            delta_x,
            delta_t,
            tau,
        }
    }
}