use lbflow::prelude::*;

fn main() {
    let n = vec![100, 100];

    let momentum_parameters = momentum::Parameters {
        n: n.clone(),
        tau: 0.9,
        delta_x: 1e-3,
        delta_t: 1.33333333e-5,
        physical_density: 998.0,
        reference_pressure: 101325.0,
        initial_density: functions::uniform_density(1.0, n.clone()),
        initial_velocity: functions::uniform_velocity(vec![0.0, 0.0], n.clone()),
        velocity_set: D2Q9,
        node_types: functions::only_fluid_nodes(n.clone()),
        boundary_conditions: vec![
            (West, momentum::bc::NoSlip),
            (East, momentum::bc::NoSlip),
            (South, momentum::bc::NoSlip),
            (
                North,
                momentum::bc::BounceBack {
                    density: 1.0,
                    velocity: vec![0.1, 0.0],
                },
            ),
        ],
        post_functions: Some(vec![momentum::PostFunction::new(
            "mean_density.csv".to_string(),
            5,
            momentum::post::compute_mean_density,
        )]),
    };

    momentum::load(momentum_parameters);
}
