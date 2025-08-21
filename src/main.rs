use lbflow::prelude::*;

fn main() {
    let n = vec![1600, 160];

    let m_params = m::Parameters {
        n: n.clone(),
        tau: 0.724637681,
        delta_x: 1e-3,
        delta_t: 1.33333333e-5,
        physical_density: 998.0,
        reference_pressure: 101325.0,
        initial_density: functions::uniform_density(1.0, n.clone()),
        initial_velocity: functions::uniform_velocity(vec![0.05, 0.0], n.clone()),
        velocity_set: D2Q9,
        node_types: functions::only_fluid_nodes(n.clone()),
        boundary_conditions: vec![
            (West, m::bc::NoSlip),
            (East, m::bc::NoSlip),
            (South, m::bc::NoSlip),
            (North, m::bc::BounceBack { density: 1.0, velocity: vec![0.1, 0.0] }),
            // (Bottom, m::bc::NoSlip),
            // (Top, m::bc::NoSlip),
        ],
        post_functions: Some(vec![m::PostFunction::new(
            "mean_density.csv".to_string(),
            5,
            m::post::compute_mean_density,
        )]),
    };

    let ps_params = ps::Parameters {
        name: "concentration".to_string(),
        tau_g: 0.72221,
        initial_concentration: functions::uniform_concentration(0.0, n.clone()),
        velocity_set: D2Q9,
        boundary_conditions: vec![
            (West, ps::bc::AntiBounceBack { concentration: 0.0 }),
            (East, ps::bc::AntiBBNoFlux),
            (South, ps::bc::AntiBounceBack { concentration: 1.0 }),
            (North, ps::bc::AntiBBNoFlux),
            // (Bottom, ps::bc::AntiBounceBack { concentration: 0.0 }),
            // (Top, ps::bc::AntiBounceBack { concentration: 0.0 }),
        ],
    };

    ps::solve(m_params, ps_params);
}
