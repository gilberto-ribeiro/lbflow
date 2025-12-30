use crate::prelude_crate::*;
use std::fmt::{self, Debug};

// -------------------------------------------------------------------------- STRUCT: Node

pub(super) struct Parameters {
    velocity_set_params: Arc<velocity_set::Parameters>,
    multiphase_params: Arc<Option<multiphase::Parameters>>,
}

impl Default for Parameters {
    fn default() -> Self {
        Parameters {
            velocity_set_params: Arc::new(velocity_set::Parameters::default()),
            multiphase_params: Arc::new(None),
        }
    }
}

impl Parameters {
    pub(super) fn new(
        velocity_set_params: Arc<velocity_set::Parameters>,
        multiphase_params: Arc<Option<multiphase::Parameters>>,
    ) -> Self {
        Parameters {
            velocity_set_params,
            multiphase_params,
        }
    }

    pub(super) fn test_default(dim: usize) -> Self {
        Parameters {
            velocity_set_params: Arc::new(velocity_set::Parameters::test_default(dim)),
            multiphase_params: Arc::new(None),
        }
    }
}

// #[derive(Debug)]
pub struct Node {
    density: RwLock<Float>,
    velocity: RwLock<Vec<Float>>,
    f: RwLock<Vec<Float>>,
    f_eq: RwLock<Vec<Float>>,
    f_star: RwLock<Vec<Float>>,
    pub(super) phi: RwLock<Option<Float>>,
    force: Arc<Option<Box<dyn Fn(&Node) -> Vec<Float> + Send + Sync>>>,
    pub(super) multiphase_parameters: Arc<Option<multiphase::Parameters>>,
    node_type: NodeType,
    index: Vec<usize>,
    coordinates: Vec<Float>,
    velocity_set_parameters: Arc<velocity_set::Parameters>,
    neighbor_nodes: RwLock<Option<HashMap<usize, Arc<Node>>>>,
    scalar_nodes: RwLock<Option<HashMap<String, Arc<passive_scalar::Node>>>>,
    bounce_back_neighbor_nodes: RwLock<Option<HashMap<usize, Arc<Node>>>>,
    bounce_back_node_status: RwLock<bool>,
    shallow_node: ShallowNode,
}

impl Debug for Node {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Node")
            .field("density", &self.get_density())
            .field("velocity", &self.get_velocity())
            .field("f", &self.get_f())
            .field("f_eq", &self.get_f_eq())
            .field("f_star", &self.get_f_star())
            .field("force", &self.get_force())
            .field("node_type", &self.node_type)
            .field("index", &self.index)
            .field("coordinates", &self.coordinates)
            .finish()
    }
}

impl Node {
    pub(super) fn new(
        density: Float,
        velocity: Vec<Float>,
        force: Arc<Option<Box<dyn Fn(&Node) -> Vec<Float> + Send + Sync>>>,
        node_type: NodeType,
        index: Vec<usize>,
        coordinates: Vec<Float>,
        parameters: Arc<Parameters>,
    ) -> Self {
        let q = parameters.velocity_set_params.q;
        Node {
            density: RwLock::new(density),
            velocity: RwLock::new(velocity.clone()),
            f: RwLock::new(vec![0.0; q]),
            f_eq: RwLock::new(vec![0.0; q]),
            f_star: RwLock::new(vec![0.0; q]),
            phi: RwLock::new(None),
            force,
            multiphase_parameters: Arc::clone(&parameters.multiphase_params),
            node_type,
            index,
            coordinates,
            velocity_set_parameters: Arc::clone(&parameters.velocity_set_params),
            neighbor_nodes: RwLock::new(None),
            scalar_nodes: RwLock::new(None),
            bounce_back_neighbor_nodes: RwLock::new(None),
            bounce_back_node_status: RwLock::new(false),
            shallow_node: ShallowNode::new(density, velocity.clone()),
        }
    }

    pub fn test_default(dim: usize) -> Self {
        match dim {
            2 => Node::new(
                1.0,
                vec![0.0, 0.0],
                Arc::new(None),
                Fluid,
                vec![3, 7],
                vec![0.035, 0.075],
                Arc::new(Parameters::test_default(2)),
            ),
            3 => Node::new(
                1.0,
                vec![0.0, 0.0, 0.0],
                Arc::new(None),
                Fluid,
                vec![3, 5, 7],
                vec![0.035, 0.055, 0.075],
                Arc::new(Parameters::test_default(3)),
            ),
            _ => panic!("Unsupported dimension: {dim}"),
        }
    }
}

impl Default for Node {
    fn default() -> Self {
        Node::new(
            1.0,
            vec![0.0, 0.0],
            Arc::new(None),
            Fluid,
            vec![0, 0],
            vec![0.0, 0.0],
            Arc::new(Parameters::default()),
        )
    }
}

impl NodeLike for Node {
    /// # Examples
    /// ```
    /// # use lbflow::momentum::Node;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_f(), vec![0.0; 9]);
    /// assert_eq!(node.get_f(), [0.0; 9]);
    /// assert_eq!(node.get_f(), &[0.0; 9]);
    /// assert_eq!(&node.get_f(), &[0.0; 9]);
    /// ```
    fn get_f(&self) -> Vec<Float> {
        self.f.read().unwrap().clone()
    }

    fn set_f(&self, f: Vec<Float>) {
        let mut f_guard = self.f.write().unwrap();
        *f_guard = f.to_vec();
    }

    fn get_f_eq(&self) -> Vec<Float> {
        self.f_eq.read().unwrap().clone()
    }

    fn set_f_eq(&self, f_eq: Vec<Float>) {
        let mut f_eq_guard = self.f_eq.write().unwrap();
        *f_eq_guard = f_eq.to_vec();
    }

    fn get_f_star(&self) -> Vec<Float> {
        self.f_star.read().unwrap().clone()
    }

    fn set_f_star(&self, f_star: Vec<Float>) {
        let mut f_star_guard = self.f_star.write().unwrap();
        *f_star_guard = f_star.to_vec();
    }

    fn get_value(&self) -> Float {
        self.get_density()
    }

    /// # Examples
    /// ## 2D
    /// ```
    /// # use lbflow::momentum::Node;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_velocity(), vec![0.0, 0.0]);
    /// assert_eq!(node.get_velocity(), [0.0, 0.0]);
    /// assert_eq!(node.get_velocity(), &[0.0, 0.0]);
    /// assert_eq!(&node.get_velocity(), &[0.0, 0.0]);
    /// ```
    /// ## 3D
    /// ```
    /// # use lbflow::momentum::Node;
    /// let node = Node::test_default(3);
    ///
    /// assert_eq!(node.get_velocity(), vec![0.0, 0.0, 0.0]);
    /// assert_eq!(node.get_velocity(), [0.0, 0.0, 0.0]);
    /// assert_eq!(node.get_velocity(), &[0.0, 0.0, 0.0]);
    /// assert_eq!(&node.get_velocity(), &[0.0, 0.0, 0.0]);
    /// ```
    fn get_velocity(&self) -> Vec<Float> {
        self.velocity.read().unwrap().clone()
    }

    fn get_vel_set_params(&self) -> &Arc<velocity_set::Parameters> {
        &self.velocity_set_parameters
    }

    /// # Examples
    /// ```
    /// # use lbflow::momentum::Node;
    /// # use lbflow::NodeType;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_node_type(), &NodeType::Fluid);
    /// assert_ne!(node.get_node_type(), &NodeType::Solid);
    /// ```
    fn get_node_type(&self) -> &NodeType {
        &self.node_type
    }

    /// # Examples
    /// ```
    /// # use lbflow::momentum::Node;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_index(), &vec![3, 7]);
    /// assert_eq!(node.get_index(), &[3, 7]);
    /// ```
    fn get_index(&self) -> &Vec<usize> {
        &self.index
    }

    fn get_coordinates(&self) -> &Vec<Float> {
        &self.coordinates
    }

    fn get_neighbor_nodes(&self) -> HashMap<usize, Arc<Node>> {
        self.neighbor_nodes
            .read()
            .unwrap()
            .as_ref()
            .cloned()
            .unwrap()
    }

    fn set_neighbor_nodes(&self, neighbor_nodes: HashMap<usize, Arc<Node>>) {
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

    fn set_bounce_back_neighbor_nodes(
        &self,
        bounce_back_neighbor_nodes: HashMap<usize, Arc<Node>>,
    ) {
        let mut bounce_back_neighbor_nodes_guard = self.bounce_back_neighbor_nodes.write().unwrap();
        *bounce_back_neighbor_nodes_guard = Some(bounce_back_neighbor_nodes);
    }

    fn is_bounce_back_node(&self) -> bool {
        *self.bounce_back_node_status.read().unwrap()
    }

    fn change_bounce_back_node_status(&self) {
        let mut status_guard = self.bounce_back_node_status.write().unwrap();
        *status_guard = !*status_guard;
    }

    fn get_scalar_nodes(&self) -> HashMap<String, Arc<passive_scalar::Node>> {
        self.scalar_nodes.read().unwrap().as_ref().cloned().unwrap()
    }

    fn compute_bgk_collision(&self, tau: Float) {
        let mut f_star = kernel::bgk_collision(
            &self.get_f(),
            &self.get_f_eq(),
            tau,
            self.get_vel_set_params(),
        );
        if let Some(force) = self.get_force().as_deref() {
            let source_term = kernel::momentum_source_term(
                &self.get_velocity(),
                force,
                tau,
                self.get_vel_set_params(),
            );
            f_star
                .iter_mut()
                .zip(source_term.iter())
                .for_each(|(f_star_i, source_term_i)| {
                    *f_star_i += *source_term_i;
                });
        };
        self.set_f_star(f_star);
    }

    fn update_shallow_node(&self) {
        let density = self.get_density();
        let velocity = self.get_velocity();
        self.get_shallow_node().set_density(density);
        self.get_shallow_node().set_velocity(velocity);
    }
}

impl Node {
    /// # Examples
    /// ```
    /// # use lbflow::momentum::Node;
    /// let node = Node::test_default(2);
    ///
    /// assert_eq!(node.get_density(), 1.0);
    /// assert_eq!(&node.get_density(), &1.0);
    /// ```
    pub fn get_density(&self) -> Float {
        *self.density.read().unwrap()
    }

    fn set_density(&self, density: Float) {
        let mut density_guard = self.density.write().unwrap();
        *density_guard = density;
    }

    fn set_velocity(&self, velocity: Vec<Float>) {
        let mut velocity_guard = self.velocity.write().unwrap();
        *velocity_guard = velocity;
    }

    fn get_force(&self) -> Option<Vec<Float>> {
        let body_force = self.force.as_ref().as_ref().map(|f| f(self));
        let shan_chen_force = self.compute_shan_chen_force();
        match (body_force, shan_chen_force) {
            (Some(bf), Some(scf)) => {
                let combined_force: Vec<Float> = bf
                    .iter()
                    .zip(scf.iter())
                    .map(|(bf_x, scf_x)| bf_x + scf_x)
                    .collect();
                Some(combined_force)
            }
            (Some(bf), None) => Some(bf),
            (None, Some(scf)) => Some(scf),
            (None, None) => None,
        }
    }

    fn get_shallow_node(&self) -> &ShallowNode {
        &self.shallow_node
    }

    pub(crate) fn append_scalar_node(&self, scalar_name: &str, node: Arc<passive_scalar::Node>) {
        self.scalar_nodes
            .write()
            .unwrap()
            .get_or_insert_with(HashMap::new)
            .insert(scalar_name.to_string(), node);
    }

    pub(super) fn compute_density(&self) {
        let f = self.get_f();
        let density = f.iter().sum::<Float>();
        self.set_density(density);
        self.compute_phi();
    }

    pub(super) fn compute_velocity(&self, explicit_computation: bool) {
        let f = self.get_f();
        let density = self.get_density();
        let vel_set_params = self.get_vel_set_params();
        let d = vel_set_params.get_d();
        let c = vel_set_params.get_c();
        let mut velocity = match (
            explicit_computation,
            vel_set_params.get_velocity_computation(),
        ) {
            (true, Some(velocity_computation)) => velocity_computation(density, f),
            (_, _) => {
                let mut velocity = Vec::with_capacity(d);
                (0..d).for_each(|x| {
                    velocity.push(
                        f.iter()
                            .zip(c.iter())
                            .map(|(f_i, c_i)| f_i * ((*c_i)[x] as Float))
                            .sum::<Float>()
                            / density,
                    );
                });
                velocity
            }
        };
        if let Some(force) = self.get_force() {
            velocity
                .iter_mut()
                .zip(force.iter())
                .for_each(|(u_x, f_x)| *u_x += 0.5 * DELTA_T * f_x / density);
        };
        self.set_velocity(velocity);
    }

    pub(super) fn compute_node_residuals(&self) -> (Float, Vec<Float>) {
        let density = self.get_density();
        let velocity = self.get_velocity();
        let shallow_density = self.get_shallow_node().get_density();
        let shallow_velocity = self.get_shallow_node().get_velocity();
        let node_residual_density = density - shallow_density;
        let node_residual_velocity = velocity
            .iter()
            .zip(shallow_velocity.iter())
            .map(|(u_x, su_x)| u_x - su_x)
            .collect::<Vec<Float>>();
        (node_residual_density, node_residual_velocity)
    }
}
// ------------------------------------------------------------------- STRUCT: ShallowNode

#[derive(Debug)]
struct ShallowNode {
    density: RwLock<Float>,
    velocity: RwLock<Vec<Float>>,
}

impl ShallowNode {
    fn new(density: Float, velocity: Vec<Float>) -> Self {
        ShallowNode {
            density: RwLock::new(density),
            velocity: RwLock::new(velocity),
        }
    }
}

impl ShallowNode {
    fn get_density(&self) -> Float {
        *self.density.read().unwrap()
    }

    fn set_density(&self, density: Float) {
        let mut density_guard = self.density.write().unwrap();
        *density_guard = density;
    }

    fn get_velocity(&self) -> Vec<Float> {
        self.velocity.read().unwrap().clone()
    }

    fn set_velocity(&self, velocity: Vec<Float>) {
        let mut velocity_guard = self.velocity.write().unwrap();
        *velocity_guard = velocity;
    }
}

#[cfg(test)]
mod tests {
    use super::super::Lattice;
    use super::*;

    #[test]
    fn test_get_f_eq_d2q9() {
        let node = Node::test_default(2);

        assert_eq!(node.get_f_eq(), vec![0.0; 9]);
        assert_eq!(node.get_f_eq(), [0.0; 9]);
        assert_eq!(node.get_f_eq(), &[0.0; 9]);
        assert_eq!(&node.get_f_eq(), &[0.0; 9]);
    }

    #[test]
    fn test_get_f_star_d2q9() {
        let node = Node::test_default(2);

        assert_eq!(node.get_f_star(), vec![0.0; 9]);
        assert_eq!(node.get_f_star(), [0.0; 9]);
        assert_eq!(node.get_f_star(), &[0.0; 9]);
        assert_eq!(&node.get_f_star(), &[0.0; 9]);
    }

    #[test]
    fn test_set_density() {
        let node = Node::test_default(2);

        assert_eq!(node.get_density(), 1.0);

        node.set_density(1.1);

        assert_eq!(node.get_density(), 1.1);
        assert_ne!(node.get_density(), 1.0);

        node.set_density(1.2);

        assert_eq!(node.get_density(), 1.2);
        assert_ne!(node.get_density(), 1.1);
    }

    #[test]
    fn test_set_velocity_d2q9() {
        let node = Node::test_default(2);

        assert_eq!(node.get_velocity(), vec![0.0, 0.0]);

        node.set_velocity(vec![0.1, 0.0]);
        let mut velocity = node.get_velocity();
        velocity[1] = 0.2;

        node.set_velocity(velocity);

        assert_eq!(node.get_velocity(), vec![0.1, 0.2]);

        let mut velocity = Vec::with_capacity(2);
        let mut u_x = 0.3;
        for _ in 0..2 {
            velocity.push(u_x);
            u_x += 0.1;
        }

        node.set_velocity(velocity);

        let actual = node.get_velocity();
        let target = [0.3, 0.4];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }

        let mut velocity = vec![0.0; 2];
        let mut u_x = 0.5;
        for elem in velocity.iter_mut() {
            *elem = u_x;
            u_x += 0.1;
        }

        node.set_velocity(velocity);

        let actual = node.get_velocity();
        let target = [0.5, 0.6];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_set_velocity_d3q19() {
        let node = Node::test_default(3);

        assert_eq!(node.get_velocity(), vec![0.0, 0.0, 0.0]);

        node.set_velocity(vec![0.1, 0.2, 0.0]);

        assert_eq!(node.get_velocity(), vec![0.1, 0.2, 0.0]);

        let mut velocity = node.get_velocity();
        velocity[2] = 0.3;

        node.set_velocity(velocity);

        assert_eq!(node.get_velocity(), vec![0.1, 0.2, 0.3]);
    }

    #[test]
    fn test_set_f_d2q9() {
        let node = Node::test_default(2);

        assert_eq!(node.get_f(), vec![0.0; 9]);

        node.set_f(vec![0.1, 0.0, 0.3, 0.0, 0.5, 0.0, 0.7, 0.0, 0.9]);

        let mut f = node.get_f();
        f[1] = 0.2;
        f[3] = 0.4;
        f[5] = 0.6;
        f[7] = 0.8;

        node.set_f(f);

        assert_eq!(
            node.get_f(),
            vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        );

        let mut f = Vec::with_capacity(9);
        let mut f_i = 0.9;
        for _ in 0..9 {
            f.push(f_i);
            f_i -= 0.1;
        }

        node.set_f(f);

        let actual = node.get_f();
        let target = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }

        let mut f = vec![0.0; 9];
        let mut f_i = 0.05;
        for elem in f.iter_mut() {
            *elem = f_i;
            f_i += 0.05;
        }

        node.set_f(f);

        let actual = node.get_f();
        let target = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_set_f_d3q19() {
        let node = Node::test_default(3);

        assert_eq!(node.get_f(), vec![0.0; 19]);

        node.set_f(vec![
            0.4, 0.3, 0.2, 0.0, 0.2, 0.3, 0.4, 0.0, 0.6, 0.7, 0.6, 0.0, 0.4, 0.3, 0.2, 0.0, 0.2,
            0.3, 0.4,
        ]);

        let mut f = node.get_f();
        f[3] = 0.1;
        f[7] = 0.5;
        f[11] = 0.5;
        f[15] = 0.1;

        node.set_f(f);

        assert_eq!(
            node.get_f(),
            vec![
                0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
                0.2, 0.3, 0.4,
            ]
        );
    }

    #[test]
    fn test_set_f_eq_d2q9() {
        let node = Node::test_default(2);

        assert_eq!(node.get_f_eq(), vec![0.0; 9]);

        node.set_f_eq(vec![0.1, 0.0, 0.3, 0.0, 0.5, 0.0, 0.7, 0.0, 0.9]);

        let mut f_eq = node.get_f_eq();
        f_eq[1] = 0.2;
        f_eq[3] = 0.4;
        f_eq[5] = 0.6;
        f_eq[7] = 0.8;

        node.set_f_eq(f_eq);

        assert_eq!(
            node.get_f_eq(),
            vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        );

        let mut f_eq = Vec::with_capacity(9);
        let mut f_eq_i = 0.9;
        for _ in 0..9 {
            f_eq.push(f_eq_i);
            f_eq_i -= 0.1;
        }

        node.set_f_eq(f_eq);

        let actual = node.get_f_eq();
        let target = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }

        let mut f_eq = vec![0.0; 9];
        let mut f_eq_i = 0.05;
        for elem in f_eq.iter_mut() {
            *elem = f_eq_i;
            f_eq_i += 0.05;
        }

        node.set_f_eq(f_eq);

        let actual = node.get_f_eq();
        let target = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_set_f_eq_d3q19() {
        let node = Node::test_default(3);

        assert_eq!(node.get_f_eq(), vec![0.0; 19]);

        node.set_f_eq(vec![
            0.4, 0.3, 0.2, 0.0, 0.2, 0.3, 0.4, 0.0, 0.6, 0.7, 0.6, 0.0, 0.4, 0.3, 0.2, 0.0, 0.2,
            0.3, 0.4,
        ]);

        let mut f_eq = node.get_f_eq();
        f_eq[3] = 0.1;
        f_eq[7] = 0.5;
        f_eq[11] = 0.5;
        f_eq[15] = 0.1;

        node.set_f_eq(f_eq);

        assert_eq!(
            node.get_f_eq(),
            vec![
                0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
                0.2, 0.3, 0.4,
            ]
        );
    }

    #[test]
    fn test_set_f_star_d2q9() {
        let node = Node::test_default(2);

        assert_eq!(node.get_f_star(), vec![0.0; 9]);

        node.set_f_star(vec![0.1, 0.0, 0.3, 0.0, 0.5, 0.0, 0.7, 0.0, 0.9]);

        let mut f_star = node.get_f_star();
        f_star[1] = 0.2;
        f_star[3] = 0.4;
        f_star[5] = 0.6;
        f_star[7] = 0.8;

        node.set_f_star(f_star);

        assert_eq!(
            node.get_f_star(),
            vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        );

        let mut f_star = Vec::with_capacity(9);
        let mut f_star_i = 0.9;
        for _ in 0..9 {
            f_star.push(f_star_i);
            f_star_i -= 0.1;
        }

        node.set_f_star(f_star);

        let actual = node.get_f_star();
        let target = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }

        let mut f_star = vec![0.0; 9];
        let mut f_star_i = 0.05;
        for elem in f_star.iter_mut() {
            *elem = f_star_i;
            f_star_i += 0.05;
        }

        node.set_f_star(f_star);

        let actual = node.get_f_star();
        let target = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_set_f_star_d3q19() {
        let node = Node::test_default(3);

        assert_eq!(node.get_f_star(), vec![0.0; 19]);

        node.set_f_star(vec![
            0.4, 0.3, 0.2, 0.0, 0.2, 0.3, 0.4, 0.0, 0.6, 0.7, 0.6, 0.0, 0.4, 0.3, 0.2, 0.0, 0.2,
            0.3, 0.4,
        ]);

        let mut f_star = node.get_f_star();
        f_star[3] = 0.1;
        f_star[7] = 0.5;
        f_star[11] = 0.5;
        f_star[15] = 0.1;

        node.set_f_star(f_star);

        assert_eq!(
            node.get_f_star(),
            vec![
                0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
                0.2, 0.3, 0.4,
            ]
        );
    }

    #[test]
    fn test_get_neighbor_node_d2q9_inner_node() {
        let momentum_lattice = Lattice::_test_default(2);
        let node = momentum_lattice._get_node_by_index(&vec![3, 7]);

        assert_eq!(node.get_neighbor_node(0).get_index(), &vec![3, 7]);
        assert_eq!(node.get_neighbor_node(1).get_index(), &vec![4, 7]);
        assert_eq!(node.get_neighbor_node(2).get_index(), &vec![3, 8]);
        assert_eq!(node.get_neighbor_node(3).get_index(), &vec![2, 7]);
        assert_eq!(node.get_neighbor_node(4).get_index(), &vec![3, 6]);
        assert_eq!(node.get_neighbor_node(5).get_index(), &vec![4, 8]);
        assert_eq!(node.get_neighbor_node(6).get_index(), &vec![2, 8]);
        assert_eq!(node.get_neighbor_node(7).get_index(), &vec![2, 6]);
        assert_eq!(node.get_neighbor_node(8).get_index(), &vec![4, 6]);
    }

    #[test]
    fn test_get_neighbor_node_d2q9_edge_node() {
        let momentum_lattice = Lattice::_test_default(2);
        let node = momentum_lattice._get_node_by_index(&vec![9, 9]);

        assert_eq!(node.get_neighbor_node(0).get_index(), &vec![9, 9]);
        assert_eq!(node.get_neighbor_node(1).get_index(), &vec![0, 9]);
        assert_eq!(node.get_neighbor_node(2).get_index(), &vec![9, 0]);
        assert_eq!(node.get_neighbor_node(3).get_index(), &vec![8, 9]);
        assert_eq!(node.get_neighbor_node(4).get_index(), &vec![9, 8]);
        assert_eq!(node.get_neighbor_node(5).get_index(), &vec![0, 0]);
        assert_eq!(node.get_neighbor_node(6).get_index(), &vec![8, 0]);
        assert_eq!(node.get_neighbor_node(7).get_index(), &vec![8, 8]);
        assert_eq!(node.get_neighbor_node(8).get_index(), &vec![0, 8]);

        assert_eq!(
            node.get_neighbor_node(3).get_neighbor_node(2).get_index(),
            &vec![8, 0]
        );
    }

    #[test]
    fn test_get_neighbor_node_d3q19_inner_node() {
        let momentum_lattice = Lattice::_test_default(3);
        let node = momentum_lattice._get_node_by_index(&vec![3, 5, 7]);

        assert_eq!(node.get_neighbor_node(0).get_index(), &vec![3, 5, 7]);
        assert_eq!(node.get_neighbor_node(1).get_index(), &vec![4, 5, 7]);
        assert_eq!(node.get_neighbor_node(7).get_index(), &vec![4, 6, 7]);
        assert_eq!(node.get_neighbor_node(9).get_index(), &vec![4, 5, 8]);
    }

    #[test]
    fn test_get_neighbor_node_d3q19_edge_node() {
        let momentum_lattice = Lattice::_test_default(3);
        let node = momentum_lattice._get_node_by_index(&vec![9, 9, 9]);

        assert_eq!(node.get_neighbor_node(0).get_index(), &vec![9, 9, 9]);
        assert_eq!(node.get_neighbor_node(1).get_index(), &vec![0, 9, 9]);
        assert_eq!(node.get_neighbor_node(7).get_index(), &vec![0, 0, 9]);
        assert_eq!(node.get_neighbor_node(9).get_index(), &vec![0, 9, 0]);
    }

    #[test]
    fn test_get_neighbor_nodes_d2q9() {
        let momentum_lattice = Lattice::_test_default(2);
        let node = momentum_lattice._get_node_by_index(&vec![3, 7]);

        let neighbor_nodes = node.get_neighbor_nodes();

        assert_eq!(neighbor_nodes[&0].get_index(), &vec![3, 7]);
        assert_eq!(neighbor_nodes[&1].get_index(), &vec![4, 7]);
        assert_eq!(neighbor_nodes[&2].get_index(), &vec![3, 8]);
        assert_eq!(neighbor_nodes[&3].get_index(), &vec![2, 7]);
        assert_eq!(neighbor_nodes[&4].get_index(), &vec![3, 6]);
        assert_eq!(neighbor_nodes[&5].get_index(), &vec![4, 8]);
        assert_eq!(neighbor_nodes[&6].get_index(), &vec![2, 8]);
        assert_eq!(neighbor_nodes[&7].get_index(), &vec![2, 6]);
        assert_eq!(neighbor_nodes[&8].get_index(), &vec![4, 6]);
    }

    #[test]
    fn test_get_neighbor_nodes_d3q19() {
        let momentum_lattice = Lattice::_test_default(3);
        let node = momentum_lattice._get_node_by_index(&vec![3, 5, 7]);

        let neighbor_nodes = node.get_neighbor_nodes();

        assert_eq!(neighbor_nodes[&0].get_index(), &vec![3, 5, 7]);
        assert_eq!(neighbor_nodes[&1].get_index(), &vec![4, 5, 7]);
        assert_eq!(neighbor_nodes[&7].get_index(), &vec![4, 6, 7]);
        assert_eq!(neighbor_nodes[&9].get_index(), &vec![4, 5, 8]);
    }

    #[test]
    fn test_get_shallow_node() {
        let node = Node::test_default(2);

        let shallow_node = node.get_shallow_node();

        assert_eq!(shallow_node.get_density(), 1.0);
        assert_eq!(shallow_node.get_velocity(), vec![0.0, 0.0]);
    }

    #[test]
    fn test_compute_density_d2q9() {
        let node = Node::test_default(2);

        node.compute_density();

        assert!((node.get_density() - 0.0).abs() < 1e-12);

        node.set_f(vec![0.1; 9]);

        node.compute_density();

        assert!((node.get_density() - 0.9).abs() < 1e-12);

        node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);

        node.compute_density();

        assert!((node.get_density() - 4.5).abs() < 1e-12);
    }

    #[test]
    fn test_compute_density_d3q19() {
        let node = Node::test_default(3);

        node.compute_density();

        assert!((node.get_density() - 0.0).abs() < 1e-12);

        node.set_f(vec![0.1; 19]);

        node.compute_density();

        assert!((node.get_density() - 1.9).abs() < 1e-12);
    }

    #[test]
    fn test_compute_velocity_d2q9_implicit_computation() {
        let node = Node::test_default(2);

        node.compute_velocity(false);

        let actual = node.get_velocity();
        let target = [0.0; 2];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }

        node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);

        node.compute_velocity(false);

        let actual = node.get_velocity();
        let target = [-0.2, -0.6];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_compute_velocity_d2q9_explicit_computation() {
        let node = Node::test_default(2);
        node.set_f(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
        node.set_density(0.5);

        node.compute_velocity(true);

        let actual = node.get_velocity();
        let target = [-0.4, -1.2];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_compute_velocity_d3q19_implicit_computation() {
        let node = Node::test_default(3);

        node.compute_velocity(false);

        let actual = node.get_velocity();
        let target = [0.0; 3];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }

        node.set_f(vec![
            0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2,
            0.3, 0.4,
        ]);

        node.compute_velocity(false);

        let actual = node.get_velocity();
        let target = [0.1, -0.3, 0.3];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_compute_velocity_d3q19_explicit_computation() {
        let node = Node::test_default(3);
        node.set_f(vec![
            0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2,
            0.3, 0.4,
        ]);
        node.set_density(0.5);

        node.compute_velocity(true);

        let actual = node.get_velocity();
        let target = [0.2, -0.6, 0.6];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_compute_bgk_collision_d2q9() {
        let node = Node::test_default(2);
        node.set_f(vec![0.1; 9]);
        node.set_f_eq(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
        let tau = 0.5;

        node.compute_bgk_collision(tau);

        let actual = node.get_f_star();
        let target = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_compute_equilibrium_d2q9() {
        let node = Node::test_default(2);
        node.set_f(vec![0.1; 9]);
        node.compute_density(); // density = 0.9

        node.compute_equilibrium();

        let actual = node.get_f_eq();
        let target = [0.4, 0.1, 0.1, 0.1, 0.1, 0.025, 0.025, 0.025, 0.025];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_compute_streaming_d2q9_case_1() {
        let momentum_lattice = Lattice::_test_default(2);
        let node = momentum_lattice._get_node_by_index(&vec![3, 7]);
        node.get_neighbor_node(0).set_f_star(vec![0.1; 9]);
        node.get_neighbor_node(1).set_f_star(vec![0.2; 9]);
        node.get_neighbor_node(2).set_f_star(vec![0.3; 9]);
        node.get_neighbor_node(3).set_f_star(vec![0.4; 9]);
        node.get_neighbor_node(4).set_f_star(vec![0.5; 9]);
        node.get_neighbor_node(5).set_f_star(vec![0.6; 9]);
        node.get_neighbor_node(6).set_f_star(vec![0.7; 9]);
        node.get_neighbor_node(7).set_f_star(vec![0.8; 9]);
        node.get_neighbor_node(8).set_f_star(vec![0.9; 9]);

        node.compute_streaming();

        assert_eq!(
            node.get_f(),
            vec![0.1, 0.4, 0.5, 0.2, 0.3, 0.8, 0.9, 0.6, 0.7]
        );
    }

    #[test]
    fn test_compute_streaming_d2q9_case_2() {
        let momentum_lattice = Lattice::_test_default(2);
        let node = momentum_lattice._get_node_by_index(&vec![3, 7]);
        node.get_neighbor_node(0)
            .set_f_star(vec![1.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
        node.get_neighbor_node(1)
            .set_f_star(vec![0.1, 0.2, 0.3, 1.0, 0.5, 0.6, 0.7, 0.8, 0.9]);
        node.get_neighbor_node(2)
            .set_f_star(vec![0.1, 0.2, 0.3, 0.4, 1.0, 0.6, 0.7, 0.8, 0.9]);
        node.get_neighbor_node(3)
            .set_f_star(vec![0.1, 1.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
        node.get_neighbor_node(4)
            .set_f_star(vec![0.1, 0.2, 1.0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);
        node.get_neighbor_node(5)
            .set_f_star(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 0.9]);
        node.get_neighbor_node(6)
            .set_f_star(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]);
        node.get_neighbor_node(7)
            .set_f_star(vec![0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 0.7, 0.8, 0.9]);
        node.get_neighbor_node(8)
            .set_f_star(vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0, 0.8, 0.9]);

        node.compute_streaming();

        assert_eq!(node.get_f(), vec![1.0; 9]);
    }

    #[test]
    fn test_compute_node_residuals() {
        let node = Node::test_default(2);
        node.set_density(0.95);
        node.set_velocity(vec![0.15, -0.25]);

        let (node_residual_density, node_residual_velocity) = node.compute_node_residuals();

        assert!((node_residual_density - (-0.05)).abs() < 1e-12);
        let actual = node_residual_velocity;
        let target = [0.15, -0.25];
        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_update_shallow_node() {
        let node = Node::test_default(2);
        let shallow_node = node.get_shallow_node();
        node.set_density(0.9);
        node.set_velocity(vec![0.1, -0.1]);

        node.update_shallow_node();

        assert_eq!(shallow_node.get_density(), 0.9);
        assert_eq!(shallow_node.get_velocity(), vec![0.1, -0.1]);
    }

    #[test]
    fn test_get_density_shallow_node() {
        let shallow_node = ShallowNode::new(1.0, vec![0.0, 0.0]);

        assert_eq!(shallow_node.get_density(), 1.0);
    }

    #[test]
    fn test_get_velocity_shallow_node() {
        let shallow_node = ShallowNode::new(1.0, vec![0.0, 0.0]);

        assert_eq!(shallow_node.get_velocity(), vec![0.0, 0.0]);
    }

    #[test]
    fn test_set_density_shallow_node() {
        let shallow_node = ShallowNode::new(1.0, vec![0.0, 0.0]);

        shallow_node.set_density(0.9);

        assert_eq!(shallow_node.get_density(), 0.9);

        shallow_node.set_density(1.1);

        assert_eq!(shallow_node.get_density(), 1.1);
    }

    #[test]
    fn test_set_velocity_shallow_node() {
        let shallow_node = ShallowNode::new(1.0, vec![0.0, 0.0]);

        shallow_node.set_velocity(vec![0.1, -0.1]);

        assert_eq!(shallow_node.get_velocity(), vec![0.1, -0.1]);

        shallow_node.set_velocity(vec![0.3, 0.7]);

        assert_eq!(shallow_node.get_velocity(), vec![0.3, 0.7]);
    }
}
