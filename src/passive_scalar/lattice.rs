use super::Node;
use super::Parameters;
use super::Residuals;
use super::bc::BoundaryCondition;
use super::bc::InnerBoundaryCondition::{self, *};
use crate::prelude_crate::*;
use crate::{FACES_2D, FACES_3D};
use rayon::prelude::*;

#[derive(Debug)]
pub struct Lattice<'a> {
    scalar_name: &'a str,
    momentum_lattice: Arc<momentum::Lattice>,
    nodes: Vec<Arc<Node>>,
    collision_operator: Arc<CollisionOperator>,
    velocity_set_parameters: Arc<velocity_set::Parameters>,
    fluid_nodes: Vec<Arc<Node>>,
    boundary_nodes: HashMap<BoundaryFace, Vec<Arc<Node>>>,
    bounce_back_nodes: Vec<Arc<Node>>,
    boundary_conditions: HashMap<BoundaryFace, BoundaryCondition>,
    inner_boundary_condition: InnerBoundaryCondition,
    residuals: RwLock<Residuals>,
}

impl<'a> Lattice<'a> {
    pub(super) fn new(params: Parameters<'a>, momentum_lattice: Arc<momentum::Lattice>) -> Self {
        let source_value = Arc::new(params.source_value);
        let collision_operator = Arc::new(params.collision_operator);
        let velocity_set = params.velocity_set;
        println!("Selecting velocity set for the lattice: {velocity_set:?}\n");
        let velocity_set_parameters = Arc::new(velocity_set.get_velocity_set_parameters());

        let d = &velocity_set_parameters.d;
        let c = &(velocity_set_parameters.c);
        let number_of_nodes = momentum_lattice.get_number_of_nodes();

        let n = momentum_lattice.get_n();
        let initial_scalar_value = &params.initial_scalar_value.generate(n, params.scalar_name);

        let adsorption_params = Arc::new(params.adsorption_parameters);
        let node_parameters = Arc::new(super::node::Parameters::new(
            Arc::clone(&velocity_set_parameters),
            adsorption_params,
        ));

        let nodes = momentum_lattice
            .get_nodes()
            .iter()
            .enumerate()
            .map(|(i, node)| {
                let initial_scalar_value = initial_scalar_value[i];
                Arc::new(Node::new(
                    initial_scalar_value,
                    Arc::clone(&source_value),
                    Arc::clone(&node_parameters),
                    Arc::clone(node),
                ))
            })
            .collect::<Vec<Arc<Node>>>();

        momentum_lattice
            .get_nodes()
            .iter()
            .zip(nodes.iter())
            .for_each(|(m_node, ps_node)| {
                m_node.append_scalar_node(params.scalar_name, Arc::clone(ps_node));
            });

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
            let index = node.get_index();
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
            let index = node.get_index();
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
            .filter(|node| matches!(node.get_node_type(), Fluid))
            .cloned()
            .collect::<Vec<Arc<Node>>>();

        let boundary_nodes = boundary_faces
            .iter()
            .zip(limit_index.iter().zip(dimension_index.iter()))
            .map(|(face, (&limit, &x))| {
                let face_nodes = fluid_nodes
                    .iter()
                    .filter(|node| node.get_index()[x] == limit)
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
                .filter(|(_, neighbor_node)| matches!(neighbor_node.get_node_type(), Solid))
                .map(|(i, neighbor_node)| (*i, Arc::clone(neighbor_node)))
                .collect::<HashMap<usize, Arc<Node>>>();
            if !bounce_back_neighbor_nodes.is_empty() {
                node.set_bounce_back_neighbor_nodes(bounce_back_neighbor_nodes);
                bounce_back_nodes.push(Arc::clone(node));
                node.change_bounce_back_node_status();
            }
        });

        Lattice {
            scalar_name: params.scalar_name,
            momentum_lattice,
            nodes,
            collision_operator: Arc::clone(&collision_operator),
            velocity_set_parameters,
            fluid_nodes,
            boundary_nodes,
            boundary_conditions: HashMap::from_iter(params.boundary_conditions),
            inner_boundary_condition: params.inner_boundary_condition,
            bounce_back_nodes,
            residuals: RwLock::new(Residuals::new(0.0)),
        }
    }
}

impl<'a> Lattice<'a> {
    pub(super) fn get_momentum_lattice(&self) -> &Arc<momentum::Lattice> {
        &self.momentum_lattice
    }

    fn _get_node(&self, i: usize) -> &Arc<Node> {
        &self.nodes[i]
    }

    pub fn get_fluid_nodes(&self) -> &Vec<Arc<Node>> {
        &self.fluid_nodes
    }

    fn get_bounce_back_nodes(&self) -> &Vec<Arc<Node>> {
        &self.bounce_back_nodes
    }

    pub(super) fn get_boundary_nodes(&self) -> &HashMap<BoundaryFace, Vec<Arc<Node>>> {
        &self.boundary_nodes
    }

    pub(super) fn get_boundary_condition(
        &self,
        boundary_face: &BoundaryFace,
    ) -> &BoundaryCondition {
        self.boundary_conditions
            .get(boundary_face)
            .expect("Boundary faces not found")
    }

    pub(super) fn get_residuals(&self) -> Residuals {
        self.residuals.read().unwrap().clone()
    }

    fn set_residuals(&self, residuals: Residuals) {
        let mut residuals_guard = self.residuals.write().unwrap();
        *residuals_guard = residuals;
    }

    pub fn get_nodes(&self) -> &Vec<Arc<Node>> {
        &self.nodes
    }

    pub(super) fn get_scalar_name(&self) -> &str {
        self.scalar_name
    }

    pub(super) fn get_d(&self) -> &usize {
        &self.velocity_set_parameters.d
    }

    fn get_collision_operator(&self) -> &Arc<CollisionOperator> {
        &self.collision_operator
    }
}

impl<'a> Lattice<'a> {
    pub(super) fn initialize_nodes(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_equilibrium();
            node.set_g(node.get_g_eq());
        });
    }

    fn update_scalar_value_step(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_scalar_value();
        });
    }

    fn equilibrium_step(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_equilibrium();
        });
    }

    fn collision_step(&self) {
        match self.get_collision_operator().as_ref() {
            BGK(tau_g) => self.get_fluid_nodes().par_iter().for_each(|node| {
                node.compute_bgk_collision(*tau_g);
            }),
            TRT(omega_plus, omega_minus) => self.get_fluid_nodes().par_iter().for_each(|node| {
                node.compute_trt_collision(*omega_plus, *omega_minus);
            }),
            MRT(relaxation_vector) => self.get_fluid_nodes().par_iter().for_each(|node| {
                node.compute_mrt_collision(relaxation_vector);
            }),
        }
    }

    fn streaming_step(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_streaming();
        });
    }

    fn inner_anti_bounce_back_step(&self) {
        self.get_bounce_back_nodes().par_iter().for_each(|node| {
            node.compute_inner_anti_bounce_back();
        });
    }

    fn inner_bounce_back_step(&self) {
        self.get_bounce_back_nodes().par_iter().for_each(|node| {
            node.compute_inner_bounce_back();
        });
    }

    fn inner_boundary_condition_step(&self) {
        match self.inner_boundary_condition {
            InnerAntiBounceBack => self.inner_anti_bounce_back_step(),
            InnerBounceBack => self.inner_bounce_back_step(),
        }
    }

    pub(super) fn compute_lattice_residuals(&self) {
        let node_residuals = self
            .get_fluid_nodes()
            .par_iter()
            .map(|node| node.compute_node_residuals())
            .collect::<Vec<Float>>();
        let residuals_scalar_value = node_residuals
            .par_iter()
            .map(|residual_scalar_value| residual_scalar_value * residual_scalar_value)
            .sum::<Float>()
            .sqrt();
        self.set_residuals(Residuals {
            scalar_value: residuals_scalar_value,
        });
    }

    pub(super) fn update_shallow_nodes(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.update_shallow_node();
        });
    }

    pub(super) fn stop_condition(&self) -> bool {
        let momentum_residuals = self.get_momentum_lattice().get_residuals();
        let converged_density = momentum_residuals.get_density() <= TOLERANCE_DENSITY;
        let converged_velocity = momentum_residuals
            .get_velocity()
            .iter()
            .all(|&u_x| u_x <= TOLERANCE_VELOCITY);
        let passive_scalar_residuals = self.get_residuals();
        let converged_scalar_value =
            passive_scalar_residuals.scalar_value <= TOLERANCE_SCALAR_VALUE;
        let converged_quantities =
            converged_density && converged_velocity && converged_scalar_value;
        let min_iterations = self.get_momentum_lattice().get_time_step() > MIN_ITER;
        let max_iterations = self.get_momentum_lattice().get_time_step()
            > self
                .get_momentum_lattice()
                .get_config()
                .get_max_iterations();
        (min_iterations && converged_quantities) || max_iterations
    }

    pub(super) fn main_steps(&self) {
        self.update_scalar_value_step();
        self.equilibrium_step();
        self.collision_step();
        self.streaming_step();
        self.inner_boundary_condition_step();
        self.boundary_conditions_step();
    }
}

pub(super) struct LatticeVec<'a> {
    passive_scalar_lattices: Vec<Lattice<'a>>,
    momentum_lattice: Arc<momentum::Lattice>,
}

impl<'a> LatticeVec<'a> {
    pub(super) fn new(
        params_vec: Vec<passive_scalar::Parameters<'a>>,
        momentum_lattice: Arc<momentum::Lattice>,
    ) -> Self {
        let passive_scalar_lattices = params_vec
            .into_iter()
            .map(|params| Lattice::new(params, Arc::clone(&momentum_lattice)))
            .collect();
        LatticeVec {
            passive_scalar_lattices,
            momentum_lattice,
        }
    }
}

impl<'a> LatticeVec<'a> {
    pub(super) fn get_momentum_lattice(&self) -> &Arc<momentum::Lattice> {
        &self.momentum_lattice
    }

    pub(super) fn get_passive_scalar_lattices(&self) -> &Vec<Lattice<'a>> {
        &self.passive_scalar_lattices
    }
}

impl<'a> LatticeVec<'a> {
    pub(super) fn initialize_nodes(&self) {
        self.passive_scalar_lattices.iter().for_each(|lattice| {
            lattice.initialize_nodes();
        });
    }

    pub(super) fn main_steps(&self) {
        self.passive_scalar_lattices.iter().for_each(|lattice| {
            lattice.main_steps();
        });
    }

    pub(super) fn compute_lattice_residuals(&self) {
        self.passive_scalar_lattices.iter().for_each(|lattice| {
            lattice.compute_lattice_residuals();
        });
    }

    pub(super) fn write_data(&self) {
        self.passive_scalar_lattices.iter().for_each(|lattice| {
            lattice.write_data();
        });
    }

    pub(super) fn update_shallow_nodes(&self) {
        self.passive_scalar_lattices.iter().for_each(|lattice| {
            lattice.update_shallow_nodes();
        });
    }

    pub(super) fn stop_condition(&self) -> bool {
        let momentum_residuals = self.get_momentum_lattice().get_residuals();
        let converged_density = momentum_residuals.get_density() <= TOLERANCE_DENSITY;
        let converged_velocity = momentum_residuals
            .get_velocity()
            .iter()
            .all(|&u_x| u_x <= TOLERANCE_VELOCITY);
        let converged_scalar_values = self
            .passive_scalar_lattices
            .iter()
            .all(|lattice| lattice.get_residuals().scalar_value <= TOLERANCE_SCALAR_VALUE);
        let converged_quantities =
            converged_density && converged_velocity && converged_scalar_values;
        let min_iterations = self.get_momentum_lattice().get_time_step() > MIN_ITER;
        let max_iterations = self.get_momentum_lattice().get_time_step()
            > self
                .get_momentum_lattice()
                .get_config()
                .get_max_iterations();
        (min_iterations && converged_quantities) || max_iterations
    }
}
