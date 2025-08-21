use super::ConversionFactor;
use super::Node;
use super::Parameters;
use super::Residuals;
use super::bc::BoundaryCondition::{self, *};
use crate::prelude::*;
use crate::{FACES_2D, FACES_3D};
use rayon::prelude::*;

const TOLERANCE_CONCENTRATION: Float = 1e-7;

#[derive(Debug)]
pub struct Lattice {
    momentum_lattice: Arc<momentum::Lattice>,
    nodes: Vec<Arc<Node>>,
    velocity_set_parameters: Arc<VelocitySetParameters>,
    _conversion_factor: Arc<ConversionFactor>,
    fluid_nodes: Vec<Arc<Node>>,
    boundary_nodes: HashMap<BoundaryFace, Vec<Arc<Node>>>,
    bounce_back_nodes: Vec<Arc<Node>>,
    boundary_conditions: HashMap<BoundaryFace, BoundaryCondition>,
    residuals: RwLock<Residuals>,
}

impl Lattice {
    pub fn new(params: Parameters, momentum_lattice: Arc<momentum::Lattice>) -> Self {
        let initial_concentration = params.initial_concentration;

        let velocity_set = params.velocity_set;
        println!("Selecting velocity set for the lattice: {velocity_set:?}\n");
        let velocity_set_parameters = Arc::new(velocity_set.get_velocity_set_parameters());

        let d = &velocity_set_parameters.d;
        let c = &(velocity_set_parameters.c);
        let number_of_nodes = momentum_lattice.get_number_of_nodes();

        let conversion_factor = Arc::new(ConversionFactor::new(
            params.tau_g,
            momentum_lattice.get_conversion_factor().delta_x,
            momentum_lattice.get_conversion_factor().delta_t,
        ));

        let nodes = momentum_lattice
            .get_nodes()
            .iter()
            .enumerate()
            .map(|(i, node)| {
                let initial_concentration = initial_concentration[i];
                Arc::new(Node::new(
                    initial_concentration,
                    Arc::clone(&velocity_set_parameters),
                    Arc::clone(&conversion_factor),
                    Arc::clone(node),
                ))
            })
            .collect::<Vec<Arc<Node>>>();

        let n = momentum_lattice.get_n();
        let nx = momentum_lattice.get_nx();
        let ny = momentum_lattice.get_ny();
        let nz = match d {
            2 => 1,
            3 => momentum_lattice.get_nz(),
            _ => panic!("Unsupported dimension: {d}"),
        };

        let mut nodes_matrix: Vec<Vec<Vec<Option<Arc<Node>>>>> = vec![vec![vec![None; nx]; ny]; nz];
        nodes.iter().for_each(|node| {
            let index = node.get_momentum_node().get_index();
            let x = index[0];
            let y = index[1];
            let z = match d {
                2 => 0,
                3 => index[2],
                _ => panic!("Unsupported dimension: {d}"),
            };
            nodes_matrix[z][y][x] = Some(Arc::clone(node));
        });

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

        println!("Setting up neighbor nodes for each node in the lattice...\n");
        nodes.iter().enumerate().for_each(|(n_i, node)| {
            crate::io::progress_bar(n_i, number_of_nodes);
            let index = node.get_momentum_node().get_index();
            let neighbor_nodes = c
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
                    let x = neighbor_index[0];
                    let y = neighbor_index[1];
                    let z = match d {
                        2 => 0,
                        3 => neighbor_index[2],
                        _ => panic!("Unsupported dimension: {d}"),
                    };
                    let neighbor_node = Arc::clone(nodes_matrix[z][y][x].as_ref().unwrap());
                    (i, neighbor_node)
                })
                .collect::<HashMap<usize, Arc<Node>>>();
            node.set_neighbor_nodes(neighbor_nodes);
        });

        let fluid_nodes = nodes
            .iter()
            .filter(|node| matches!(node.get_momentum_node().get_node_type(), Fluid))
            .cloned()
            .collect::<Vec<Arc<Node>>>();

        let boundary_nodes = boundary_faces
            .iter()
            .zip(limit_index.iter().zip(dimension_index.iter()))
            .map(|(face, (&limit, &x))| {
                let face_nodes = fluid_nodes
                    .iter()
                    .filter(|node| node.get_momentum_node().get_index()[x] == limit)
                    .cloned()
                    .collect::<Vec<Arc<Node>>>();
                (*face, face_nodes)
            })
            .collect::<HashMap<BoundaryFace, Vec<Arc<Node>>>>();

        let mut bounce_back_nodes = Vec::new();
        fluid_nodes.iter().for_each(|node| {
            let bounce_back_neighbor_nodes = node
                .get_neighbor_nodes()
                .iter()
                .filter(|(_, neighbor_node)| {
                    matches!(neighbor_node.get_momentum_node().get_node_type(), Solid)
                })
                .map(|(i, neighbor_node)| (*i, Arc::clone(neighbor_node)))
                .collect::<HashMap<usize, Arc<Node>>>();
            if !bounce_back_neighbor_nodes.is_empty() {
                // node.bounce_back_neighbor_nodes
                //     .replace(Some(bounce_back_neighbor_nodes));
                node.set_bounce_back_neighbor_nodes(bounce_back_neighbor_nodes);
                bounce_back_nodes.push(Arc::clone(node));
            }
        });

        Lattice {
            momentum_lattice,
            nodes,
            velocity_set_parameters,
            _conversion_factor: Arc::clone(&conversion_factor),
            fluid_nodes,
            boundary_nodes,
            boundary_conditions: HashMap::from_iter(params.boundary_conditions),
            bounce_back_nodes,
            residuals: RwLock::new(Residuals::new(0.0)),
        }
    }
}

impl Lattice {
    pub fn get_momentum_lattice(&self) -> &Arc<momentum::Lattice> {
        &self.momentum_lattice
    }

    pub fn get_node(&self, i: usize) -> &Arc<Node> {
        &self.nodes[i]
    }

    pub fn get_fluid_nodes(&self) -> &Vec<Arc<Node>> {
        &self.fluid_nodes
    }

    pub fn get_bounce_back_nodes(&self) -> &Vec<Arc<Node>> {
        &self.bounce_back_nodes
    }

    pub fn get_boundary_nodes(&self) -> &HashMap<BoundaryFace, Vec<Arc<Node>>> {
        &self.boundary_nodes
    }

    pub fn get_boundary_condition(&self, boundary_face: &BoundaryFace) -> &BoundaryCondition {
        self.boundary_conditions
            .get(boundary_face)
            .expect("Boundary faces not found")
    }

    pub fn get_residuals(&self) -> Residuals {
        // self.residuals.borrow().clone()
        self.residuals.read().unwrap().clone()
    }

    pub fn set_residuals(&self, residuals: Residuals) {
        // self.residuals.replace(residuals);
        let mut residuals_guard = self.residuals.write().unwrap();
        *residuals_guard = residuals;
    }

    pub fn get_nodes(&self) -> &Vec<Arc<Node>> {
        &self.nodes
    }
}

impl Lattice {
    pub fn get_d(&self) -> &usize {
        &self.velocity_set_parameters.d
    }
}

impl Lattice {
    pub fn initialize_nodes(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_equilibrium();
            node.set_g(node.get_g_eq());
        });
    }

    pub fn update_concentration_step(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_concentration();
        });
    }

    pub fn equilibrium_step(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_equilibrium();
        });
    }

    pub fn bgk_collision_step(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_bgk_collision();
        });
    }

    pub fn streaming_step(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_streaming();
        });
    }

    pub fn inner_anti_bounce_back_step(&self) {
        self.get_bounce_back_nodes().par_iter().for_each(|node| {
            node.compute_inner_anti_bounce_back();
        });
    }

    pub fn boundary_conditions_step(&self) {
        self.get_boundary_nodes()
            .iter()
            .for_each(|(boundary_face, nodes)| {
                let boundary_condition = self.get_boundary_condition(boundary_face);
                let momentum_boundary_condition = self
                    .get_momentum_lattice()
                    .get_boundary_condition(boundary_face);
                let dim = *self.get_d();
                let velocity = match momentum_boundary_condition {
                    momentum::bc::NoSlip => Some(vec![0.0; dim]),
                    momentum::bc::BounceBack { velocity, .. } => Some(velocity.to_vec()),
                    momentum::bc::AntiBounceBack { .. } => None,
                    momentum::bc::Periodic => None,
                };
                match boundary_condition {
                    AntiBounceBack { concentration } => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_anti_bounce_back_bc(
                                boundary_face,
                                concentration,
                                velocity.as_deref(),
                            );
                        });
                    }
                    AntiBBNoFlux => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_no_flux_bc(boundary_face, velocity.as_deref());
                        });
                    }
                    Periodic => {}
                }
            });
    }

    pub fn compute_lattice_residuals(&self) {
        let node_residuals = self
            .get_fluid_nodes()
            .par_iter()
            .map(|node| node.compute_node_residuals())
            .collect::<Vec<Float>>();
        let residuals_concentration = node_residuals
            .par_iter()
            .map(|residual_concentration| residual_concentration * residual_concentration)
            .sum::<Float>()
            .sqrt();
        self.set_residuals(Residuals {
            concentration: residuals_concentration,
        });
    }

    pub fn update_shallow_nodes(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.update_shallow_node();
        });
    }

    pub fn stop_condition(&self) -> bool {
        let momentum_residuals = self.get_momentum_lattice().get_residuals();
        let converged_density = momentum_residuals.density <= momentum::lattice::TOLERANCE_DENSITY;
        let converged_velocity = momentum_residuals
            .velocity
            .iter()
            .all(|&u_x| u_x <= momentum::lattice::TOLERANCE_VELOCITY);
        let passive_scalar_residuals = self.get_residuals();
        let converged_concentration =
            passive_scalar_residuals.concentration <= TOLERANCE_CONCENTRATION;
        let converged_quantities =
            converged_density && converged_velocity && converged_concentration;
        let min_iterations =
            self.get_momentum_lattice().get_time_step() > momentum::lattice::MIN_ITER;
        let max_iterations = self.get_momentum_lattice().get_time_step()
            > self
                .get_momentum_lattice()
                .get_config()
                .get_max_iterations();
        (min_iterations && converged_quantities) || max_iterations
    }
}
