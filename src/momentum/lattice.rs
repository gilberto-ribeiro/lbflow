use super::bc::BoundaryCondition::{self, *};
use super::post::PostFunction;
use super::{ConversionFactor, Node, Parameters, Residuals};
use crate::prelude::*;
use crate::{FACES_2D, FACES_3D};
use rayon::prelude::*;

pub const MIN_ITER: usize = 10;
pub const TOLERANCE_DENSITY: Float = 1e-7;
pub const TOLERANCE_VELOCITY: Float = 1e-7;

// ----------------------------------------------------------------------- STRUCT: Lattice

#[derive(Debug)]
pub struct Lattice {
    nodes: Vec<Arc<Node>>,
    n: Vec<usize>,
    collision_operator: Arc<CollisionOperator>,
    velocity_set_parameters: Arc<VelocitySetParameters>,
    conversion_factor: Arc<ConversionFactor>,
    fluid_nodes: Vec<Arc<Node>>,
    boundary_nodes: HashMap<BoundaryFace, Vec<Arc<Node>>>,
    bounce_back_nodes: Vec<Arc<Node>>,
    boundary_conditions: HashMap<BoundaryFace, BoundaryCondition>,
    residuals: RwLock<Residuals>,
    time_step: RwLock<usize>,
    post_functions: Option<Vec<PostFunction>>,
    config: Config,
}

impl Lattice {
    pub fn new(config: Config, params: Parameters) -> Self {
        let initial_density = &params.initial_density;
        let initial_velocity = &params.initial_velocity;

        let collision_operator = Arc::new(params.collision_operator);
        let velocity_set = params.velocity_set;
        println!("Selecting collision operator for the lattice: {collision_operator:?}\n");
        println!("Selecting velocity set for the lattice: {velocity_set:?}\n");
        let velocity_set_parameters = Arc::new(velocity_set.get_velocity_set_parameters());

        let d = &velocity_set_parameters.d;
        let c = &(velocity_set_parameters.c);
        let n = &params.n;
        let node_types = &params.node_types;
        let physical_delta_x = params.delta_x;
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

        let conversion_factor = Arc::new(ConversionFactor::new(
            Arc::clone(&collision_operator),
            params.velocity_set,
            params.delta_x,
            params.delta_t,
            params.physical_density,
            params.reference_pressure,
        ));

        println!("Creating lattice with dimensions: {n:?}\n");
        let nodes = (0..num_nodes)
            .map(|i| {
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
                let node_type = node_types[i];
                let density = initial_density[i];
                let velocity = initial_velocity[i].clone();
                Arc::new(Node::new(
                    density,
                    velocity,
                    node_type,
                    index,
                    coordinates,
                    Arc::clone(&velocity_set_parameters),
                ))
            })
            .collect::<Vec<Arc<Node>>>();

        let nx = n[0];
        let ny = n[1];
        let nz = match d {
            2 => 1,
            3 => n[2],
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

        println!("Setting up neighbor nodes for each node in the lattice...\n");
        nodes.iter().enumerate().for_each(|(n_i, node)| {
            crate::io::progress_bar(n_i, num_nodes);
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
                    let neighbor_node = Arc::clone(nodes_matrix[z][y][x].as_ref().unwrap());
                    (i, neighbor_node)
                })
                .collect::<HashMap<usize, Arc<Node>>>();
            node.set_neighbor_nodes(neighbor_nodes);
        });

        println!("Identifying fluid nodes...\n");
        let fluid_nodes = nodes
            .iter()
            .filter(|node| matches!(node.get_node_type(), Fluid))
            .cloned()
            .collect::<Vec<Arc<Node>>>();

        println!("Identifying boundary nodes...\n");
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

        println!("Identifying bounce-back nodes...\n");
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
            }
        });

        Lattice {
            nodes,
            n: n.to_vec(),
            collision_operator: Arc::clone(&collision_operator),
            velocity_set_parameters: Arc::clone(&velocity_set_parameters),
            conversion_factor,
            fluid_nodes,
            boundary_nodes,
            boundary_conditions: HashMap::from_iter(params.boundary_conditions),
            bounce_back_nodes,
            residuals: RwLock::new(Residuals::new(0.0, vec![0.0; *d])),
            time_step: RwLock::new(0),
            post_functions: params.post_functions,
            config,
        }
    }

    pub fn test_default(dim: usize) -> Self {
        let parameters = Parameters::test_default(dim);
        let config = Config::default();
        Lattice::new(config, parameters)
    }
}

impl Default for Lattice {
    fn default() -> Self {
        let parameters = Parameters::default();
        let config = Config::default();
        Lattice::new(config, parameters)
    }
}

impl Lattice {
    /// # Examples
    /// ```
    /// # use std::sync::Arc;
    /// # use lbflow::momentum::{Lattice, Node};
    /// # use lbflow::momentum::Parameters;
    /// let momentum_lattice = Lattice::test_default(2);
    ///
    /// assert_eq!(momentum_lattice.get_node_by_index(&vec![3, 7]).get_index(), &vec![3, 7]);
    pub fn get_node_by_index(&self, index: &Vec<usize>) -> &Arc<Node> {
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
    pub fn get_boundary_nodes_by_face(&self, boundary_face: &BoundaryFace) -> &Vec<Arc<Node>> {
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

    pub fn get_nodes(&self) -> &Vec<Arc<Node>> {
        &self.nodes
    }

    /// Examples
    /// 2D
    /// ```
    /// # use std::sync::Arc;
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

    pub fn get_time_step(&self) -> usize {
        // *self.time_step.borrow()
        *self.time_step.read().unwrap()
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

    pub fn get_nx(&self) -> usize {
        self.n[0]
    }

    pub fn get_ny(&self) -> usize {
        self.n[1]
    }

    pub fn get_nz(&self) -> usize {
        match self.get_d() {
            2 => 1,
            3 => self.n[2],
            _ => panic!("Unsupported dimension: {}", self.get_d()),
        }
    }

    pub fn get_number_of_nodes(&self) -> usize {
        self.n.par_iter().product()
    }

    pub fn get_n(&self) -> &Vec<usize> {
        &self.n
    }

    pub fn get_config(&self) -> &Config {
        &self.config
    }

    pub fn get_post_functions(&self) -> &Option<Vec<PostFunction>> {
        &self.post_functions
    }
}

impl Lattice {
    pub fn get_d(&self) -> &usize {
        &self.velocity_set_parameters.d
    }

    pub fn get_conversion_factor(&self) -> &Arc<ConversionFactor> {
        &self.conversion_factor
    }

    pub fn get_collision_operator(&self) -> &Arc<CollisionOperator> {
        &self.collision_operator
    }
}

impl Lattice {
    pub fn initialize_nodes(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_equilibrium();
            node.set_f(node.get_f_eq());
        });
    }

    pub fn update_density_and_velocity_step(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_density();
            node.compute_velocity(false);
        });
    }

    pub fn equilibrium_step(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_equilibrium();
        });
    }

    pub fn collision_step(&self) {
        match self.get_collision_operator().as_ref() {
            BGK(tau) => self.get_fluid_nodes().par_iter().for_each(|node| {
                node.compute_bgk_collision(*tau);
            }),
            TRT(omega_plus, omega_minus) => self.get_fluid_nodes().par_iter().for_each(|node| {
                node.compute_trt_collision(*omega_plus, *omega_minus);
            }),
            MRT(relaxation_vector) => self.get_fluid_nodes().par_iter().for_each(|node| {
                node.compute_mrt_collision(relaxation_vector);
            }),
        }
    }

    pub fn streaming_step(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_streaming();
        });
    }

    pub fn inner_bounce_back_step(&self) {
        self.get_bounce_back_nodes().par_iter().for_each(|node| {
            node.compute_inner_bounce_back();
        });
    }

    pub fn boundary_conditions_step(&self) {
        self.get_boundary_nodes()
            .iter()
            .for_each(|(boundary_face, nodes)| {
                let boundary_condition = self.get_boundary_condition(boundary_face);
                match boundary_condition {
                    NoSlip => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_no_slip_bc(boundary_face);
                        });
                    }
                    BounceBack { density, velocity } => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_bounce_back_bc(boundary_face, density, velocity);
                        });
                    }
                    AntiBounceBack { density } => {
                        nodes.par_iter().for_each(|node| {
                            node.compute_anti_bounce_back_bc(boundary_face, density);
                        });
                    }
                    Periodic => {}
                }
            });
    }

    pub fn compute_lattice_residuals(&self) {
        let d = self.get_d();
        let node_residuals = self
            .get_fluid_nodes()
            .par_iter()
            .map(|node| node.compute_node_residuals())
            .collect::<Vec<(Float, Vec<Float>)>>();
        let residuals_density = node_residuals
            .par_iter()
            .map(|(residual_density, _)| residual_density * residual_density)
            .sum::<Float>()
            .sqrt();
        let residuals_velocity = (0..*d)
            .map(|x| {
                node_residuals
                    .par_iter()
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
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.update_shallow_node();
        });
    }

    pub fn next_time_step(&self) {
        // *self.time_step.borrow_mut() += 1;
        let mut time_step_guard = self.time_step.write().unwrap();
        *time_step_guard += 1;
    }

    pub fn stop_condition(&self) -> bool {
        let residuals = self.get_residuals();
        let converged_density = residuals.density <= TOLERANCE_DENSITY;
        let converged_velocity = residuals
            .velocity
            .iter()
            .all(|&u_x| u_x <= TOLERANCE_VELOCITY);
        let converged_quantities = converged_density && converged_velocity;
        let min_iterations = self.get_time_step() > MIN_ITER;
        let max_iterations = self.get_time_step() > self.get_config().get_max_iterations();
        (min_iterations && converged_quantities) || max_iterations
    }

    pub fn compute_post_processing(&self) {
        if let Some(post_functions) = self.get_post_functions() {
            post_functions.iter().for_each(|post_function| {
            self.write_post_processing(post_function).unwrap_or_else(|e| {
                eprintln! {"Error while writing the post-processing file {}: {}", post_function.file_name, e};
                std::process::exit(1);
            });
        });
        }
    }

    pub fn main_steps(&self) {
        if !self.get_config().freeze_momentum {
            self.update_density_and_velocity_step();
            self.equilibrium_step();
            self.collision_step();
            self.streaming_step();
            self.inner_bounce_back_step();
            self.boundary_conditions_step();
        }
    }
}
