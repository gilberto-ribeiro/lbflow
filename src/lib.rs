mod momentum;
mod thermal;

type Float = f64;

use std::rc::Rc;

pub fn run() {
    let mut momentum_lattice = Rc::new(momentum::Lattice::new(4, 4, 1.0, 0.1, 0.6));
    // momentum_lattice.nodes[0].density = 5.0;
    let thermal_lattice = Rc::new(thermal::Lattice::new(
        Rc::clone(&momentum_lattice),
        0.5,
    ));

    println!("Momentum lattice created with {} nodes.", thermal_lattice.momentum_lattice.nodes.len());
    println!("Thermal lattice created with {} thermal nodes.", thermal_lattice.nodes.len());

    println!("First thermal node: {:#?}", thermal_lattice.nodes[0]);
    println!();
    println!("First momentum node: {:#?}", momentum_lattice.nodes[0]);
}