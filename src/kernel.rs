use crate::prelude::*;

pub fn equilibrium(
    dens_conc: Float,
    velocity: &[Float],
    vel_set_params: &VelocitySetParameters,
) -> Vec<Float> {
    let q = vel_set_params.get_q();
    let c = vel_set_params.get_c();
    let w = vel_set_params.get_w();
    let mut f_eq = Vec::with_capacity(q);
    let u_dot_u = velocity.iter().map(|u_x| u_x * u_x).sum::<Float>();
    (0..q).for_each(|i| {
        let u_dot_c = velocity
            .iter()
            .zip(c[i].iter())
            .map(|(u_x, c_x)| u_x * (*c_x as Float))
            .sum::<Float>();
        f_eq.push(
            w[i] * dens_conc
                * (1.0 + u_dot_c * CS_2_INV + 0.5 * u_dot_c * u_dot_c * CS_4_INV
                    - 0.5 * u_dot_u * CS_2_INV),
        );
    });
    f_eq
}

pub fn bgk_collision(
    velocity: &[Float],
    f: &[Float],
    f_eq: &[Float],
    tau: Float,
    force: Option<&[Float]>,
    vel_set_params: &VelocitySetParameters,
) -> Vec<Float> {
    let q = vel_set_params.get_q();
    let omega = DELTA_T / tau;
    let omega_prime = 1.0 - omega;
    let mut f_star = Vec::with_capacity(q);
    (0..q).for_each(|i| {
        f_star.push(omega_prime * f[i] + omega * f_eq[i]);
    });
    if let Some(force) = force {
        let source_term = kernel::source_term(velocity, force, tau, vel_set_params);
        f_star
            .iter_mut()
            .zip(source_term.iter())
            .for_each(|(f_star_i, source_term_i)| {
                *f_star_i += *source_term_i;
            });
    };
    f_star
}

pub fn trt_collision(
    f: &[Float],
    f_eq: &[Float],
    omega_plus: Float,
    omega_minus: Float,
    vel_set_params: &VelocitySetParameters,
) -> Vec<Float> {
    let q = vel_set_params.get_q();
    let mut f_star = Vec::with_capacity(q);
    (0..q).for_each(|i| {
        let i_bar = vel_set_params.get_opposite_direction(i);
        let f_plus = 0.5 * (f[i] + f[i_bar]);
        let f_eq_plus = 0.5 * (f_eq[i] + f_eq[i_bar]);
        let f_minus = 0.5 * (f[i] - f[i_bar]);
        let f_eq_minus = 0.5 * (f_eq[i] - f_eq[i_bar]);
        f_star.push(
            f[i] - omega_plus * DELTA_T * (f_plus - f_eq_plus)
                - omega_minus * DELTA_T * (f_minus - f_eq_minus),
        );
    });
    f_star
}

pub fn mrt_collision(
    dens_conc: Float,
    velocity: &[Float],
    f: &[Float],
    f_eq: &[Float],
    relaxation_vector: &[Float],
    vel_set_params: &VelocitySetParameters,
) -> Vec<Float> {
    let q = vel_set_params.get_q();
    let mrt_matrix = vel_set_params.get_mrt_matrix();
    let mrt_inverse_matrix = vel_set_params.get_mrt_inverse_matrix();
    let mut m = Vec::with_capacity(q);
    let mut m_eq = Vec::with_capacity(q);
    let mut m_star = Vec::with_capacity(q);
    let mut f_star = Vec::with_capacity(q);
    (0..q).for_each(|k| {
        m.push(
            mrt_matrix[k]
                .iter()
                .zip(f.iter())
                .map(|(matrix_ki, f_i)| matrix_ki * f_i)
                .sum::<Float>(),
        )
    });
    match vel_set_params.get_mrt_equilibrium_moments_computation() {
        Some(&mrt_equilibrium_moments_computation) => {
            m_eq = mrt_equilibrium_moments_computation(dens_conc, velocity.to_vec());
        }
        None => {
            (0..q).for_each(|k| {
                m_eq.push(
                    mrt_matrix[k]
                        .iter()
                        .zip(f_eq.iter())
                        .map(|(matrix_ki, f_eq_i)| matrix_ki * f_eq_i)
                        .sum::<Float>(),
                )
            });
        }
    }
    (0..q).for_each(|k| {
        m_star.push(m[k] - relaxation_vector[k] * DELTA_T * (m[k] - m_eq[k]));
    });
    (0..q).for_each(|i| {
        f_star.push(
            mrt_inverse_matrix[i]
                .iter()
                .zip(m_star.iter())
                .map(|(inverse_matrix_ik, m_star_k)| inverse_matrix_ik * m_star_k)
                .sum::<Float>(),
        );
    });
    f_star
}

pub fn source_term(
    velocity: &[Float],
    force: &[Float],
    tau: Float,
    vel_set_params: &VelocitySetParameters,
) -> Vec<Float> {
    let q = vel_set_params.get_q();
    let c = vel_set_params.get_c();
    let w = vel_set_params.get_w();
    let coeff = 1.0 - 0.5 * DELTA_T / tau;
    let mut source_term = Vec::with_capacity(q);
    (0..q).for_each(|i| {
        let u_dot_c = velocity
            .iter()
            .zip(c[i].iter())
            .map(|(u_x, c_x)| u_x * (*c_x as Float))
            .sum::<Float>();
        source_term.push(
            coeff
                * w[i]
                * velocity
                    .iter()
                    .zip(force.iter().zip(c[i].iter()))
                    .map(|(u_x, (f_x, c_ix))| {
                        (CS_2_INV * (*c_ix as Float) - CS_2_INV * u_x
                            + CS_4_INV * u_dot_c * (*c_ix as Float))
                            * f_x
                    })
                    .sum::<Float>(),
        );
    });
    source_term
}
