use super::adsorption;
use crate::prelude_crate::*;
use std::fmt::{self, Debug};

pub(super) struct Parameters {
    velocity_set_params: Arc<velocity_set::Parameters>,
    adsorption_params: Arc<Option<adsorption::Parameters>>,
}

impl Default for Parameters {
    fn default() -> Self {
        Parameters {
            velocity_set_params: Arc::new(velocity_set::Parameters::default()),
            adsorption_params: Arc::new(None),
        }
    }
}

impl Parameters {
    pub(super) fn new(
        velocity_set_params: Arc<velocity_set::Parameters>,
        adsorption_params: Arc<Option<adsorption::Parameters>>,
    ) -> Self {
        Parameters {
            velocity_set_params,
            adsorption_params,
        }
    }

    pub(super) fn _test_default(dim: usize) -> Self {
        Parameters {
            velocity_set_params: Arc::new(velocity_set::Parameters::test_default(dim)),
            adsorption_params: Arc::new(None),
        }
    }
}

// #[derive(Debug)]
pub struct Node {
    momentum_node: Arc<momentum::Node>,
    scalar_value: RwLock<Float>,
    g: RwLock<Vec<Float>>,
    g_eq: RwLock<Vec<Float>>,
    g_star: RwLock<Vec<Float>>,
    pub(super) q_ads: RwLock<Option<Float>>,
    pub(super) adsorption_parameters: Arc<Option<adsorption::Parameters>>,
    source_value: Arc<Option<Box<dyn Fn(&Node) -> Float + Send + Sync>>>,
    velocity_set_parameters: Arc<velocity_set::Parameters>,
    neighbor_nodes: RwLock<Option<HashMap<usize, Arc<Node>>>>,
    bounce_back_neighbor_nodes: RwLock<Option<HashMap<usize, Arc<Node>>>>,
    bounce_back_node_status: RwLock<bool>,
    shallow_node: ShallowNode,
}

impl Debug for Node {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Node")
            .field("momentum_node", &self.momentum_node)
            .field("scalar_value", &self.get_scalar_value())
            .field("g", &self.get_g())
            .field("g_eq", &self.get_g_eq())
            .field("g_star", &self.get_g_star())
            .field("velocity_set_parameters", &self.velocity_set_parameters)
            .finish()
    }
}

impl Node {
    pub(super) fn new(
        scalar_value: Float,
        source_value: Arc<Option<Box<dyn Fn(&Node) -> Float + Send + Sync>>>,
        parameters: Arc<Parameters>,
        momentum_node: Arc<momentum::Node>,
    ) -> Self {
        let q = parameters.velocity_set_params.q;
        let q_ads_value = if let Some(ads_params) = parameters.adsorption_params.as_ref() {
            if momentum_node.is_bounce_back_node() {
                Some(ads_params.get_initial_q_ads())
            } else {
                Some(0.0)
            }
        } else {
            None
        };
        Node {
            momentum_node,
            scalar_value: RwLock::new(scalar_value),
            g: RwLock::new(vec![0.0; q]),
            g_eq: RwLock::new(vec![0.0; q]),
            g_star: RwLock::new(vec![0.0; q]),
            q_ads: RwLock::new(q_ads_value),
            source_value,
            velocity_set_parameters: Arc::clone(&parameters.velocity_set_params),
            adsorption_parameters: Arc::clone(&parameters.adsorption_params),
            neighbor_nodes: RwLock::new(None),
            bounce_back_neighbor_nodes: RwLock::new(None),
            bounce_back_node_status: RwLock::new(false),
            shallow_node: ShallowNode::new(scalar_value),
        }
    }
}

impl Node {
    pub fn get_momentum_node(&self) -> &Arc<momentum::Node> {
        &self.momentum_node
    }

    pub fn get_scalar_value(&self) -> Float {
        *self.scalar_value.read().unwrap()
    }

    fn set_scalar_value(&self, scalar_value: Float) {
        let mut scalar_value_guard = self.scalar_value.write().unwrap();
        *scalar_value_guard = scalar_value;
    }

    pub fn get_g(&self) -> Vec<Float> {
        self.g.read().unwrap().clone()
    }

    pub(super) fn set_g(&self, g: Vec<Float>) {
        let mut g_guard = self.g.write().unwrap();
        *g_guard = g;
    }

    pub(super) fn get_g_eq(&self) -> Vec<Float> {
        self.g_eq.read().unwrap().clone()
    }

    pub(super) fn set_g_eq(&self, g_eq: Vec<Float>) {
        let mut g_eq_guard = self.g_eq.write().unwrap();
        *g_eq_guard = g_eq;
    }

    pub(super) fn get_g_star(&self) -> Vec<Float> {
        self.g_star.read().unwrap().clone()
    }

    pub(super) fn set_g_star(&self, g_star: Vec<Float>) {
        let mut g_star_guard = self.g_star.write().unwrap();
        *g_star_guard = g_star;
    }

    pub fn is_bounce_back_node(&self) -> bool {
        *self.bounce_back_node_status.read().unwrap()
    }

    pub(super) fn change_bounce_back_node_status(&self) {
        let mut status_guard = self.bounce_back_node_status.write().unwrap();
        *status_guard = !*status_guard;
    }

    pub(super) fn get_neighbor_nodes(&self) -> HashMap<usize, Arc<Node>> {
        self.neighbor_nodes
            .read()
            .unwrap()
            .as_ref()
            .cloned()
            .unwrap()
    }

    pub(super) fn set_neighbor_nodes(&self, neighbor_nodes: HashMap<usize, Arc<Node>>) {
        let mut neighbor_nodes_guard = self.neighbor_nodes.write().unwrap();
        *neighbor_nodes_guard = Some(neighbor_nodes);
    }

    fn get_bounce_back_neighbor_nodes(&self) -> HashMap<usize, Arc<Node>> {
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
        let mut bounce_back_neighbor_nodes_guard = self.bounce_back_neighbor_nodes.write().unwrap();
        *bounce_back_neighbor_nodes_guard = Some(bounce_back_neighbor_nodes);
    }

    pub(super) fn get_neighbor_node(&self, i: usize) -> Arc<Node> {
        self.get_neighbor_nodes()
            .get(&i)
            .cloned()
            .expect("Neighbor node not found")
    }

    pub fn get_scalar_node(&self, scalar_name: String) -> Arc<passive_scalar::Node> {
        self.get_momentum_node().get_scalar_node(scalar_name)
    }

    fn get_shallow_node(&self) -> &ShallowNode {
        &self.shallow_node
    }

    fn get_source_value(&self) -> Option<Float> {
        let original_source_value = self.source_value.as_ref().as_ref().map(|f| f(self));
        let adsorption_source_value = self.compute_adsorption_source_value();
        match (original_source_value, adsorption_source_value) {
            (Some(v1), Some(v2)) => Some(v1 + v2),
            (Some(v), None) | (None, Some(v)) => Some(v),
            (None, None) => None,
        }
    }
}

impl Node {
    pub(super) fn compute_scalar_value(&self) {
        let g = self.get_g();
        let mut scalar_value = g.iter().sum::<Float>();
        if let Some(source_value) = self.get_source_value() {
            scalar_value += 0.5 * DELTA_T * source_value;
        };
        self.set_scalar_value(scalar_value);
    }

    pub(super) fn compute_equilibrium(&self) {
        let g_eq = kernel::equilibrium(
            self.get_scalar_value(),
            &self.get_momentum_node().get_velocity(),
            self.get_velocity_set_parameters(),
        );
        self.set_g_eq(g_eq);
    }

    pub(super) fn compute_bgk_collision(&self, tau_g: Float) {
        let mut g_star = kernel::bgk_collision(
            &self.get_g(),
            &self.get_g_eq(),
            tau_g,
            self.get_velocity_set_parameters(),
        );
        if let Some(source_value) = self.get_source_value() {
            let source_term = kernel::passive_scalar_source_term(
                source_value,
                tau_g,
                self.get_velocity_set_parameters(),
            );
            g_star
                .iter_mut()
                .zip(source_term.iter())
                .for_each(|(g_star_i, source_term_i)| {
                    *g_star_i += *source_term_i;
                });
        }
        self.set_g_star(g_star);
    }

    pub(super) fn compute_trt_collision(&self, omega_plus: Float, omega_minus: Float) {
        let g_star = kernel::trt_collision(
            &self.get_g(),
            &self.get_g_eq(),
            omega_plus,
            omega_minus,
            self.get_velocity_set_parameters(),
        );
        self.set_g_star(g_star);
    }

    pub(super) fn compute_mrt_collision(&self, relaxation_vector: &[Float]) {
        let g_star = kernel::mrt_collision(
            self.get_scalar_value(),
            &self.get_momentum_node().get_velocity(),
            &self.get_g(),
            &self.get_g_eq(),
            relaxation_vector,
            self.get_velocity_set_parameters(),
        );
        self.set_g_star(g_star);
    }

    pub(super) fn compute_streaming(&self) {
        let vel_set_params = self.get_velocity_set_parameters();
        let q = vel_set_params.get_q();
        let mut g = vec![0.0; q];
        self.get_neighbor_nodes()
            .iter()
            .for_each(|(i, neighbor_node)| {
                let i_bar = vel_set_params.get_opposite_direction(*i);
                g[i_bar] = neighbor_node.get_g_star()[i_bar];
            });
        self.set_g(g);
    }

    pub(super) fn compute_inner_anti_bounce_back(&self) {
        let mut g = self.get_g();
        let g_star = self.get_g_star();
        let vel_set_params = self.get_velocity_set_parameters();
        let w = vel_set_params.get_w();
        self.get_bounce_back_neighbor_nodes()
            .iter()
            .for_each(|(&i, node)| {
                let i_bar = vel_set_params.get_opposite_direction(i);
                let wall_scalar_value = node.get_scalar_value();
                g[i_bar] = -g_star[i] + 2.0 * w[i] * wall_scalar_value;
            });
        self.set_g(g);
    }

    pub(super) fn compute_inner_bounce_back(&self) {
        let mut g = self.get_g();
        let g_star = self.get_g_star();
        let vel_set_params = self.get_velocity_set_parameters();
        self.get_bounce_back_neighbor_nodes()
            .iter()
            .for_each(|(&i, _)| {
                let i_bar = vel_set_params.get_opposite_direction(i);
                g[i_bar] = g_star[i];
            });
        self.set_g(g);
    }

    pub(super) fn compute_node_residuals(&self) -> Float {
        let scalar_value = self.get_scalar_value();
        let shallow_scalar_value = self.get_shallow_node().get_scalar_value();
        scalar_value - shallow_scalar_value
    }

    pub(super) fn update_shallow_node(&self) {
        let scalar_value = self.get_scalar_value();
        self.get_shallow_node().set_scalar_value(scalar_value);
    }

    pub(super) fn get_velocity_set_parameters(&self) -> &Arc<velocity_set::Parameters> {
        &self.velocity_set_parameters
    }
}

#[derive(Debug)]
struct ShallowNode {
    scalar_value: RwLock<Float>,
}

impl ShallowNode {
    fn new(scalar_value: Float) -> Self {
        ShallowNode {
            scalar_value: RwLock::new(scalar_value),
        }
    }
}

impl ShallowNode {
    fn get_scalar_value(&self) -> Float {
        *self.scalar_value.read().unwrap()
    }

    fn set_scalar_value(&self, scalar_value: Float) {
        let mut scalar_guard = self.scalar_value.write().unwrap();
        *scalar_guard = scalar_value;
    }
}
