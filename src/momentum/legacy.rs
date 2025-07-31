use super::{Lattice, Node, Residuals, ShallowNode};
use crate::constants::*;
use crate::velocity_set::d2q9::{C, D, Q, W};

impl Node {
    fn legacy_compute_density(&self) {
        let f = self.get_f();
        let density = f.iter().sum::<Float>();
        self.set_density(density);
    }

    fn legacy_compute_velocity(&self) {
        let density = self.get_density();
        let f = self.get_f();
        let mut velocity = self.get_velocity();
        velocity[0] = (1.0 / density) * (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]);
        velocity[1] = (1.0 / density) * (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]);
        self.set_velocity(velocity);
    }

    fn legacy_compute_moments(&self) {
        self.legacy_compute_density();
        self.legacy_compute_velocity();
    }

    fn legacy_compute_equilibrium(&self, explicit_version: bool) {
        let density = self.get_density();
        let velocity = self.get_velocity();
        let mut f_eq = self.get_f_eq();
        let u_x = velocity[0];
        let u_y = velocity[1];
        let u_x_2 = u_x * u_x;
        let u_y_2 = u_y * u_y;
        let u_2 = u_x_2 + u_y_2;
        if explicit_version {
            let ux_uy = u_x * u_y;
            let mut coeff = (2.0 * density) / 9.0;
            f_eq[0] = coeff * (2.0 - 3.0 * u_2);
            coeff = density / 18.0;
            f_eq[1] = coeff * (2.0 + 6.0 * u_x + 9.0 * u_x_2 - 3.0 * u_2);
            f_eq[2] = coeff * (2.0 + 6.0 * u_y + 9.0 * u_y_2 - 3.0 * u_2);
            f_eq[3] = coeff * (2.0 - 6.0 * u_x + 9.0 * u_x_2 - 3.0 * u_2);
            f_eq[4] = coeff * (2.0 - 6.0 * u_y + 9.0 * u_y_2 - 3.0 * u_2);
            coeff = density / 36.0;
            f_eq[5] = coeff * (1.0 + 3.0 * (u_x + u_y) + 9.0 * ux_uy + 3.0 * u_2);
            f_eq[6] = coeff * (1.0 - 3.0 * (u_x - u_y) - 9.0 * ux_uy + 3.0 * u_2);
            f_eq[7] = coeff * (1.0 - 3.0 * (u_x + u_y) + 9.0 * ux_uy + 3.0 * u_2);
            f_eq[8] = coeff * (1.0 + 3.0 * (u_x - u_y) - 9.0 * ux_uy + 3.0 * u_2);
        } else {
            for i in 0..Q {
                let cx = C[i][0] as Float;
                let cy = C[i][1] as Float;
                let u_dot_c = u_x * cx + u_y * cy;
                f_eq[i] = W[i]
                    * density
                    * (1.0 + CS_2_INV * u_dot_c + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                        - 0.5 * CS_2_INV * u_2);
            }
        }
        self.set_f_eq(f_eq);
    }

    fn legacy_compute_bgk_collision(&self, omega: Float, omega_prime: Float) {
        let f = self.get_f();
        let f_eq = self.get_f_eq();
        let mut f_star = self.get_f_star();
        for i in 0..Q {
            f_star[i] = omega_prime * f[i] + omega * f_eq[i];
        }
        self.set_f_star(f_star);
    }
}

impl Lattice {
    fn legacy_streaming_step(&self) {
        for x in 0..self.n[0] {
            for y in 0..self.n[1] {
                for i in 0..Q {
                    let [cx, cy] = C[i];
                    let new_x = ((x as i32) + cx).rem_euclid(self.n[0] as i32) as usize;
                    let new_y = ((y as i32) + cy).rem_euclid(self.n[1] as i32) as usize;
                    let node = self.get_node_by_index(&vec![x, y]);
                    let neighbor_node = self.get_node_by_index(&vec![new_x, new_y]);
                    let f_star = node.get_f_star();
                    let mut f = neighbor_node.get_f();
                    f[i] = f_star[i];
                    neighbor_node.set_f(f);
                }
            }
        }
    }

    fn legacy_compute_lattice_residuals(&self) -> Residuals {
        let old_lattice = self
            .nodes
            .iter()
            .map(|node| node.get_shallow_node().clone())
            .collect::<Vec<ShallowNode>>();
        let density = self
            .nodes
            .iter()
            .zip(old_lattice.iter())
            .map(|(node, old_node)| (node.get_density() - old_node.get_density()).powi(2))
            .sum::<Float>()
            .sqrt();
        let velocity_x = self
            .nodes
            .iter()
            .zip(old_lattice.iter())
            .map(|(node, old_node)| (node.get_velocity()[0] - old_node.get_velocity()[0]).powi(2))
            .sum::<Float>()
            .sqrt();
        let velocity_y = self
            .nodes
            .iter()
            .zip(old_lattice.iter())
            .map(|(node, old_node)| (node.get_velocity()[1] - old_node.get_velocity()[1]).powi(2))
            .sum::<Float>()
            .sqrt();
        Residuals {
            density,
            velocity: vec![velocity_x, velocity_y],
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::constants::DELTA_T;
    use crate::momentum::Node;
    use rand::Rng;

    #[test]
    fn compare_legacy_compute_bgk_collision() {
        let legacy_node = Node::test_default(2);

        let node = Node::test_default(2);

        let random_tau: f64 = rand::rng().random_range(0.5..1.0);

        let omega = DELTA_T / random_tau;
        let omega_prime = 1.0 - omega;

        let random_f: Vec<f64> = (0..9).map(|_| rand::rng().random_range(0.0..1.0)).collect();
        let random_f_eq: Vec<f64> = (0..9).map(|_| rand::rng().random_range(0.0..1.0)).collect();

        legacy_node.set_f(random_f.clone());
        legacy_node.set_f_eq(random_f_eq.clone());

        node.set_f(random_f.clone());
        node.set_f_eq(random_f_eq.clone());

        legacy_node.legacy_compute_bgk_collision(omega, omega_prime);
        node.compute_bgk_collision(omega, omega_prime);
        let actual = legacy_node.get_f_star();
        let target = node.get_f_star();

        for (a, b) in actual.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }
}
