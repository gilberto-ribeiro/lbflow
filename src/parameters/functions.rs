use crate::NodeType;

pub fn fluid_and_solid_nodes_by_bounce_back_map() -> Vec<NodeType> {
    crate::io::read_bounce_back_map().expect("Failed to read bounce-back map")
}

pub fn only_fluid_nodes(n: Vec<usize>) -> Vec<NodeType> {
    let num_nodes = n.iter().product();
    vec![NodeType::Fluid; num_nodes]
}
