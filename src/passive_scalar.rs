use crate::prelude::*;
use crate::{FACES_2D, FACES_3D};
use rayon::prelude::*;

#[derive(Debug)]
pub struct Node {
    momentum_node: Arc<momentum::Node>,
    concentration: RwLock<Float>,
    g: RwLock<Vec<Float>>,
    g_eq: RwLock<Vec<Float>>,
    g_star: RwLock<Vec<Float>>,
    velocity_set_parameters: Arc<VelocitySetParameters>,
    neighbor_nodes: RwLock<Option<HashMap<usize, Arc<Node>>>>,
    bounce_back_neighbor_nodes: RwLock<Option<HashMap<usize, Arc<Node>>>>,
}

pub struct Parameters {
    pub tau_g: Float,
    pub initial_concentration: Float,
    pub velocity_set: VelocitySet,
}

#[derive(Debug)]
pub struct Lattice {
    momentum_lattice: Arc<momentum::Lattice>,
    nodes: Vec<Arc<Node>>,
    tau_g: Float,
    velocity_set_parameters: Arc<VelocitySetParameters>,
    fluid_nodes: Vec<Arc<Node>>,
    boundary_nodes: HashMap<BoundaryFace, Vec<Arc<Node>>>,
    bounce_back_nodes: Vec<Arc<Node>>,
}

impl Node {
    pub fn new(
        concentration: Float,
        velocity_set_parameters: Arc<VelocitySetParameters>,
        momentum_node: Arc<momentum::Node>,
    ) -> Self {
        let q = velocity_set_parameters.q;
        Node {
            momentum_node,
            concentration: RwLock::new(concentration),
            g: RwLock::new(vec![0.0; q]),
            g_eq: RwLock::new(vec![0.0; q]),
            g_star: RwLock::new(vec![0.0; q]),
            velocity_set_parameters,
            neighbor_nodes: RwLock::new(None),
            bounce_back_neighbor_nodes: RwLock::new(None),
        }
    }
}

impl Node {
    pub fn get_momentum_node(&self) -> &Arc<momentum::Node> {
        &self.momentum_node
    }

    pub fn get_density(&self) -> Float {
        self.get_momentum_node().get_density()
    }

    pub fn get_velocity(&self) -> Vec<Float> {
        self.get_momentum_node().get_velocity()
    }

    pub fn get_concentration(&self) -> Float {
        // self.concentration.borrow().clone()
        *self.concentration.read().unwrap()
    }

    pub fn set_concentration(&self, concentration: Float) {
        // self.concentration.replace(concentration);
        let mut concentration_lock = self.concentration.write().unwrap();
        *concentration_lock = concentration;
    }

    pub fn get_g(&self) -> Vec<Float> {
        // self.g.borrow().clone()
        self.g.read().unwrap().clone()
    }

    fn set_g(&self, g: Vec<Float>) {
        // self.g.replace(g);
        let mut g_lock = self.g.write().unwrap();
        *g_lock = g;
    }

    pub fn get_g_eq(&self) -> Vec<Float> {
        // self.g_eq.borrow().clone()
        self.g_eq.read().unwrap().clone()
    }

    fn set_g_eq(&self, g_eq: Vec<Float>) {
        // self.g_eq.replace(g_eq);
        let mut g_eq_lock = self.g_eq.write().unwrap();
        *g_eq_lock = g_eq;
    }

    pub fn get_g_star(&self) -> Vec<Float> {
        // self.g_star.borrow().clone()
        self.g_star.read().unwrap().clone()
    }

    fn set_g_star(&self, g_star: Vec<Float>) {
        // self.g_star.replace(g_star);
        let mut g_star_lock = self.g_star.write().unwrap();
        *g_star_lock = g_star;
    }

    pub fn get_index(&self) -> &Vec<usize> {
        &self.get_momentum_node().get_index()
    }

    pub fn get_coordinates(&self) -> &Vec<Float> {
        &self.get_momentum_node().get_coordinates()
    }

    pub fn get_node_type(&self) -> &NodeType {
        &self.get_momentum_node().get_node_type()
    }

    pub fn get_neighbor_nodes(&self) -> HashMap<usize, Arc<Node>> {
        // self.neighbor_nodes.borrow().as_ref().cloned().unwrap()
        self.neighbor_nodes
            .read()
            .unwrap()
            .as_ref()
            .cloned()
            .unwrap()
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
}

impl Node {
    pub fn compute_concentration(&self) {
        let g = self.get_g();
        let concentration = g.iter().sum::<Float>();
        self.set_concentration(concentration);
    }

    pub fn compute_equilibrium(&self) {
        let concentration = self.get_concentration();
        let velocity = self.get_velocity();
        let q = self.get_q();
        let c = self.get_c();
        let w = self.get_w();
        let mut g_eq = Vec::with_capacity(*q);
        (0..*q).for_each(|i| {
            let u_dot_c = velocity
                .iter()
                .zip(c[i].iter())
                .map(|(u_x, c_x)| u_x * (*c_x as Float))
                .sum::<Float>();
            let u_dot_u = velocity.iter().map(|u_x| u_x * u_x).sum::<Float>();
            g_eq.push(
                w[i] * concentration
                    * (1.0 + u_dot_c * CS_2_INV + 0.5 * u_dot_c * u_dot_c * CS_4_INV
                        - 0.5 * u_dot_u * CS_2_INV),
            );
        });
        self.set_g_eq(g_eq);
    }

    pub fn compute_bgk_collision(&self, omega_g: Float, omega_g_prime: Float) {
        let g = self.get_g();
        let g_eq = self.get_g_eq();
        let q = self.get_q();
        let mut f_star = Vec::with_capacity(*q);
        (0..*q).for_each(|i| {
            f_star.push(omega_g_prime * g[i] + omega_g * g_eq[i]);
        });
        self.set_g_star(f_star);
    }

    pub fn compute_streaming(&self) {
        let q = self.get_q();
        let mut g = vec![0.0; *q];
        self.get_neighbor_nodes()
            .iter()
            .for_each(|(i, neighbor_node)| {
                let i_bar = self.get_opposite_direction(*i);
                g[i_bar] = neighbor_node.get_g_star()[i_bar];
            });
        self.set_g(g);
    }

    pub fn compute_inner_anti_bounce_back(&self, wall_concentration: Float) {
        let mut g = self.get_g();
        let g_star = self.get_g_star();
        let w = self.get_w();
        self.get_bounce_back_neighbor_nodes()
            .iter()
            .for_each(|(&i, _)| {
                let i_bar = self.get_opposite_direction(i);
                g[i_bar] = -g_star[i] + 2.0 * w[i] * wall_concentration;
            });
        self.set_g(g);
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

    pub fn get_velocity_computation(&self) -> Option<fn(Float, Vec<Float>) -> Vec<Float>> {
        self.get_velocity_set_parameters().velocity_computation
    }

    pub fn get_opposite_direction(&self, direction: usize) -> usize {
        self.get_q_bar()[direction]
    }
}

impl Lattice {
    pub fn new(parameters: Parameters, momentum_lattice: Arc<momentum::Lattice>) -> Self {
        let initial_concentration = parameters.initial_concentration;
        let velocity_set = parameters.velocity_set;
        println!("Selecting velocity set for the lattice: {velocity_set:?}\n");
        let velocity_set_parameters = Arc::new(velocity_set.get_velocity_set_parameters());

        let d = &velocity_set_parameters.d;
        let c = &(velocity_set_parameters.c);
        let number_of_nodes = momentum_lattice.get_number_of_nodes();

        let nodes: Vec<Arc<Node>> = momentum_lattice
            .get_nodes()
            .iter()
            .map(|node| {
                Arc::new(Node::new(
                    initial_concentration,
                    Arc::clone(&velocity_set_parameters),
                    Arc::clone(node),
                ))
            })
            .collect();

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
            let mut neighbor_nodes_guard = node.neighbor_nodes.write().unwrap();
            *neighbor_nodes_guard = Some(neighbors);
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
                // node.bounce_back_neighbor_nodes
                //     .replace(Some(bounce_back_neighbor_nodes));
                let mut bounce_back_neighbor_nodes_guard =
                    node.bounce_back_neighbor_nodes.write().unwrap();
                *bounce_back_neighbor_nodes_guard = Some(bounce_back_neighbor_nodes);
                bounce_back_nodes.push(Arc::clone(node));
            }
        });

        Lattice {
            momentum_lattice,
            nodes,
            tau_g: parameters.tau_g,
            velocity_set_parameters,
            fluid_nodes,
            boundary_nodes,
            bounce_back_nodes,
        }
    }
}

impl Lattice {
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
}

impl Lattice {
    pub fn initialize_nodes(&self) {
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_equilibrium();
            node.set_g(node.get_g_eq());
        });
    }

    pub fn update_concentration(&self) {
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
        let omega_g = DELTA_T / self.tau_g;
        let omega_g_prime = 1.0 - omega_g;
        self.get_fluid_nodes().par_iter().for_each(|node| {
            node.compute_bgk_collision(omega_g, omega_g_prime);
        });
    }

    // pub fn streaming_step(&self) {
    //     self.get_fluid_nodes().par_iter().for_each(|node| {
    //         node.compute_streaming();
    //     });
    // }

    // pub fn inner_bounce_back_step(&self) {
    //     self.get_bounce_back_nodes().par_iter().for_each(|node| {
    //         node.compute_inner_bounce_back();
    //     });
    // }

    // pub fn boundary_conditions_step(&self) {
    //     self.get_boundary_nodes()
    //         .iter()
    //         .for_each(|(boundary_face, nodes)| {
    //             let boundary_condition = self.get_boundary_condition(boundary_face);
    //             match boundary_condition {
    //                 NoSlip => {
    //                     nodes.par_iter().for_each(|node| {
    //                         node.compute_no_slip_bc(boundary_face);
    //                     });
    //                 }
    //                 BounceBack { density, velocity } => {
    //                     nodes.par_iter().for_each(|node| {
    //                         node.compute_bounce_back_bc(boundary_face, density, velocity);
    //                     });
    //                 }
    //                 AntiBounceBack { density } => {
    //                     nodes.par_iter().for_each(|node| {
    //                         node.compute_anti_bounce_back_bc(boundary_face, density);
    //                     });
    //                 }
    //                 Periodic => {}
    //             }
    //         });
    // }

    // pub fn compute_lattice_residuals(&self) {
    //     let d = self.get_d();
    //     let node_residuals = self
    //         .get_fluid_nodes()
    //         .par_iter()
    //         .map(|node| node.compute_node_residuals())
    //         .collect::<Vec<(Float, Vec<Float>)>>();
    //     let residuals_density = node_residuals
    //         .par_iter()
    //         .map(|(residual_density, _)| residual_density * residual_density)
    //         .sum::<Float>()
    //         .sqrt();
    //     let residuals_velocity = (0..*d)
    //         .map(|x| {
    //             node_residuals
    //                 .par_iter()
    //                 .map(|(_, residual_velocity)| residual_velocity[x] * residual_velocity[x])
    //                 .sum::<Float>()
    //                 .sqrt()
    //         })
    //         .collect::<Vec<Float>>();
    //     self.set_residuals(Residuals {
    //         density: residuals_density,
    //         velocity: residuals_velocity,
    //     });
    // }

    // pub fn update_shallow_nodes(&self) {
    //     self.get_fluid_nodes().par_iter().for_each(|node| {
    //         node.update_shallow_node();
    //     });
    // }
}
