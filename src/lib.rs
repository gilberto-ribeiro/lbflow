// ------------------------------------------------------------------------------- MODULES

mod cli;
mod constants;
mod io;
mod kernel;
pub mod momentum;
pub mod passive_scalar;
pub mod prelude;
mod prelude_crate;
mod velocity_set;

// ------------------------------------------------------------------------------- IMPORTS

use prelude_crate::*;

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum BoundaryFace {
    West = 0,
    East = 1,
    South = 2,
    North = 3,
    Bottom = 4,
    Top = 5,
}

const FACES_2D: [BoundaryFace; 4] = [West, East, South, North];

const FACES_3D: [BoundaryFace; 6] = [West, East, South, North, Bottom, Top];

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum NodeType {
    Fluid = 0,
    Solid = 1,
}

#[derive(Debug)]
pub enum CollisionOperator {
    BGK(Float),
    TRT(Float, Float),
    MRT(Vec<Float>),
}

impl Default for CollisionOperator {
    fn default() -> Self {
        BGK(0.5)
    }
}

pub(crate) trait NodeLike {
    fn get_f(&self) -> Vec<Float>;

    fn set_f(&self, f: Vec<Float>);

    fn get_f_eq(&self) -> Vec<Float>;

    fn set_f_eq(&self, f_eq: Vec<Float>);

    fn get_f_star(&self) -> Vec<Float>;

    fn set_f_star(&self, f_star: Vec<Float>);

    fn get_value(&self) -> Float;

    // fn set_value(&self, value: Float);

    fn get_velocity(&self) -> Vec<Float>;

    // fn set_velocity(&self, velocity: Vec<Float>);

    fn get_vel_set_params(&self) -> &Arc<velocity_set::Parameters>;

    fn get_node_type(&self) -> &NodeType;

    fn get_index(&self) -> &Vec<usize>;

    fn get_coordinates(&self) -> &Vec<Float>;

    fn get_neighbor_nodes(&self) -> HashMap<usize, Arc<Self>>;

    fn get_neighbor_node(&self, i: usize) -> Arc<Self> {
        self.get_neighbor_nodes()
            .get(&i)
            .cloned()
            .expect("Neighbor node not found")
    }

    fn set_neighbor_nodes(&self, neighbor_nodes: HashMap<usize, Arc<Self>>);

    fn get_bounce_back_neighbor_nodes(&self) -> HashMap<usize, Arc<Self>>;

    fn set_bounce_back_neighbor_nodes(&self, bounce_back_neighbor_nodes: HashMap<usize, Arc<Self>>);

    fn is_bounce_back_node(&self) -> bool;

    fn change_bounce_back_node_status(&self);

    fn get_scalar_nodes(&self) -> HashMap<String, Arc<passive_scalar::Node>>;

    fn get_scalar_node(&self, scalar_name: String) -> Arc<passive_scalar::Node> {
        self.get_scalar_nodes()
            .get(&scalar_name)
            .cloned()
            .expect("Scalar node not found")
    }

    fn compute_equilibrium(&self) {
        let f_eq = kernel::equilibrium(
            self.get_value(),
            &self.get_velocity(),
            self.get_vel_set_params(),
        );
        self.set_f_eq(f_eq);
    }

    fn compute_bgk_collision(&self, tau: Float);

    fn compute_trt_collision(&self, omega_plus: Float, omega_minus: Float) {
        let f_star = kernel::trt_collision(
            &self.get_f(),
            &self.get_f_eq(),
            omega_plus,
            omega_minus,
            self.get_vel_set_params(),
        );
        self.set_f_star(f_star);
    }

    fn compute_mrt_collision(&self, relaxation_vector: &[Float]) {
        let f_star = kernel::mrt_collision(
            self.get_value(),
            &self.get_velocity(),
            &self.get_f(),
            &self.get_f_eq(),
            relaxation_vector,
            self.get_vel_set_params(),
        );
        self.set_f_star(f_star);
    }

    fn compute_streaming(&self) {
        let vel_set_params = self.get_vel_set_params();
        let q = vel_set_params.get_q();
        let mut f = vec![0.0; q];
        self.get_neighbor_nodes()
            .iter()
            .for_each(|(i, neighbor_node)| {
                let i_bar = vel_set_params.get_opposite_direction(*i);
                f[i_bar] = neighbor_node.get_f_star()[i_bar];
            });
        self.set_f(f);
    }

    fn compute_inner_bounce_back(&self) {
        let mut f = self.get_f();
        let f_star = self.get_f_star();
        let vel_set_params = self.get_vel_set_params();
        self.get_bounce_back_neighbor_nodes()
            .iter()
            .for_each(|(&i, _)| {
                let i_bar = vel_set_params.get_opposite_direction(i);
                f[i_bar] = f_star[i];
            });
        self.set_f(f);
    }

    fn update_shallow_node(&self);
}

pub fn solve(
    momentum_params: momentum::Parameters,
    passive_scalar_params: Option<Vec<passive_scalar::Parameters>>,
) {
    match passive_scalar_params {
        None => {
            momentum::solve(momentum_params);
        }
        Some(passive_scalar_params_vec) => match passive_scalar_params_vec.len() {
            0 => {
                panic!("Must provide at least one passive_scalar_parameters inside the vector.");
            }
            // 1 => {
            //     let passive_scalar_params = &passive_scalar_params_vec[0];
            //     passive_scalar::solve(momentum_params, passive_scalar_params);
            // }
            _ => {
                passive_scalar::solve_vec(momentum_params, passive_scalar_params_vec);
            }
        },
    }
}
