use super::ConversionFactor;
use crate::prelude::*;

#[derive(Debug)]
pub struct Node {
    momentum_node: Arc<momentum::Node>,
    concentration: RwLock<Float>,
    g: RwLock<Vec<Float>>,
    g_eq: RwLock<Vec<Float>>,
    g_star: RwLock<Vec<Float>>,
    velocity_set_parameters: Arc<VelocitySetParameters>,
    conversion_factor: Arc<ConversionFactor>,
    neighbor_nodes: RwLock<Option<HashMap<usize, Arc<Node>>>>,
    bounce_back_neighbor_nodes: RwLock<Option<HashMap<usize, Arc<Node>>>>,
    shallow_node: ShallowNode,
}

impl Node {
    pub fn new(
        concentration: Float,
        velocity_set_parameters: Arc<VelocitySetParameters>,
        conversion_factor: Arc<ConversionFactor>,
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
            conversion_factor,
            neighbor_nodes: RwLock::new(None),
            bounce_back_neighbor_nodes: RwLock::new(None),
            shallow_node: ShallowNode::new(concentration),
        }
    }
}

impl Node {
    pub fn get_momentum_node(&self) -> &Arc<momentum::Node> {
        &self.momentum_node
    }

    pub fn get_concentration(&self) -> Float {
        *self.concentration.read().unwrap()
    }

    pub(super) fn set_concentration(&self, concentration: Float) {
        let mut concentration_lock = self.concentration.write().unwrap();
        *concentration_lock = concentration;
    }

    pub fn get_g(&self) -> Vec<Float> {
        self.g.read().unwrap().clone()
    }

    pub(super) fn set_g(&self, g: Vec<Float>) {
        let mut g_lock = self.g.write().unwrap();
        *g_lock = g;
    }

    pub fn get_g_eq(&self) -> Vec<Float> {
        self.g_eq.read().unwrap().clone()
    }

    pub(super) fn set_g_eq(&self, g_eq: Vec<Float>) {
        let mut g_eq_lock = self.g_eq.write().unwrap();
        *g_eq_lock = g_eq;
    }

    pub fn get_g_star(&self) -> Vec<Float> {
        self.g_star.read().unwrap().clone()
    }

    pub(super) fn set_g_star(&self, g_star: Vec<Float>) {
        let mut g_star_lock = self.g_star.write().unwrap();
        *g_star_lock = g_star;
    }

    pub fn get_neighbor_nodes(&self) -> HashMap<usize, Arc<Node>> {
        self.neighbor_nodes
            .read()
            .unwrap()
            .as_ref()
            .cloned()
            .unwrap()
    }

    pub(super) fn set_neighbor_nodes(&self, neighbor_nodes: HashMap<usize, Arc<Node>>) {
        // self.neighbor_nodes.replace(Some(neighbors));
        let mut neighbor_nodes_guard = self.neighbor_nodes.write().unwrap();
        *neighbor_nodes_guard = Some(neighbor_nodes);
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

    pub(super) fn set_bounce_back_neighbor_nodes(
        &self,
        bounce_back_neighbor_nodes: HashMap<usize, Arc<Node>>,
    ) {
        // node.bounce_back_neighbor_nodes
        //     .replace(Some(bounce_back_neighbor_nodes));
        let mut bounce_back_neighbor_nodes_guard = self.bounce_back_neighbor_nodes.write().unwrap();
        *bounce_back_neighbor_nodes_guard = Some(bounce_back_neighbor_nodes);
    }

    pub fn get_neighbor_node(&self, i: usize) -> Arc<Node> {
        self.get_neighbor_nodes()
            .get(&i)
            .cloned()
            .expect("Neighbor node not found")
    }

    pub fn get_shallow_node(&self) -> &ShallowNode {
        &self.shallow_node
    }
}

impl Node {
    pub fn get_conversion_factor(&self) -> &Arc<ConversionFactor> {
        &self.conversion_factor
    }

    pub fn get_tau_g(&self) -> Float {
        self.get_conversion_factor().tau_g
    }
}

impl Node {
    pub fn compute_concentration(&self) {
        let g = self.get_g();
        let concentration = g.iter().sum::<Float>();
        self.set_concentration(concentration);
    }

    pub fn compute_equilibrium(&self) {
        let g_eq = kernel::equilibrium(
            self.get_concentration(),
            &self.get_momentum_node().get_velocity(),
            self.get_velocity_set_parameters(),
        );
        self.set_g_eq(g_eq);
    }

    pub fn compute_bgk_collision(&self, tau_g: Float) {
        let g_star = kernel::bgk_collision(
            &self.get_momentum_node().get_velocity(),
            &self.get_g(),
            &self.get_g_eq(),
            tau_g,
            None,
            self.get_velocity_set_parameters(),
        );
        self.set_g_star(g_star);
    }

    pub fn compute_trt_collision(&self, omega_plus: Float, omega_minus: Float) {
        let g_star = kernel::trt_collision(
            &self.get_g(),
            &self.get_g_eq(),
            omega_plus,
            omega_minus,
            self.get_velocity_set_parameters(),
        );
        self.set_g_star(g_star);
    }

    pub fn compute_mrt_collision(&self, relaxation_vector: &[Float]) {
        let g_star = kernel::mrt_collision(
            self.get_concentration(),
            &self.get_momentum_node().get_velocity(),
            &self.get_g(),
            &self.get_g_eq(),
            relaxation_vector,
            self.get_velocity_set_parameters(),
        );
        self.set_g_star(g_star);
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

    pub fn compute_inner_anti_bounce_back(&self) {
        let mut g = self.get_g();
        let g_star = self.get_g_star();
        let w = self.get_w();
        self.get_bounce_back_neighbor_nodes()
            .iter()
            .for_each(|(&i, node)| {
                let i_bar = self.get_opposite_direction(i);
                let wall_concentration = node.get_concentration();
                g[i_bar] = -g_star[i] + 2.0 * w[i] * wall_concentration;
            });
        self.set_g(g);
    }

    pub fn compute_node_residuals(&self) -> Float {
        let concentration = self.get_concentration();
        let shallow_concentration = self.get_shallow_node().get_concentration();
        concentration - shallow_concentration
    }

    pub fn update_shallow_node(&self) {
        let concentration = self.get_concentration();
        self.get_shallow_node().set_concentration(concentration);
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

    pub fn get_opposite_direction(&self, direction: usize) -> usize {
        self.get_q_bar()[direction]
    }
}

#[derive(Debug)]
pub struct ShallowNode {
    pub concentration: RwLock<Float>,
}

impl ShallowNode {
    pub fn new(density: Float) -> Self {
        ShallowNode {
            concentration: RwLock::new(density),
        }
    }
}

impl ShallowNode {
    pub fn get_concentration(&self) -> Float {
        *self.concentration.read().unwrap()
    }

    pub fn set_concentration(&self, concentration: Float) {
        let mut concentration_guard = self.concentration.write().unwrap();
        *concentration_guard = concentration;
    }
}
