use crate::constants::G;
use crate::orbit::simulate_system;
use crate::system::System;
use ndarray::{array, concatenate, Array1, Array2, Axis};

/// Apply impulse to a system
///
/// The system initial position and velocities as well as the gravitational parameter are changed according to the effect of the impulse. The effect of impulse are a mass ejection affecting mu and a velocity change affecting the new initial velocity to be set.
///
/// **Inputs**:
///
/// system: the system to which the impulse is applied
///
/// rr: the position at which the impulse is applied
///
/// vv: the velocity of the body being simulated at the time the impulse is applied
///
/// delta_vv: change in the velocity from the impulse
///
/// delta_m: change in the mass of the simulated body from the impulse
///
/// **Output**: a System struct with new initial position and velocity and a new gravitational parameter. Note that this function does not change the starting or ending times attached of the system.
pub fn apply_impulse_to_system(
    system: System,
    rr: Array1<f64>,
    vv: Array1<f64>,
    delta_vv: Array1<f64>,
    delta_m: f64,
) -> System {
    let mut new_system: System = system.clone();
    new_system
        .f64_parameters
        .insert("mu".to_string(), system.f64_parameters["mu"] + delta_m / G);
    new_system
        .array_parameters
        .insert("vv0".to_string(), vv + delta_vv);
    new_system.array_parameters.insert("rr0".to_string(), rr);
    new_system
}

///Helper to create patch system from last simulation to point to impulse time point
fn create_patch_system_orbit_to_impulse(
    system: System,
    impulse_t: f64,
    simulated_orbit: &Array2<f64>,
) -> System {
    let sim_len: usize = simulated_orbit.dim().1 as usize;
    let mut new_system: System = system.clone();
    new_system.start_time = simulated_orbit[[0, sim_len - 1]];
    new_system.end_time = impulse_t;
    new_system.array_parameters.insert(
        "rr0".to_string(),
        array!(
            simulated_orbit[[1, sim_len - 1]],
            simulated_orbit[[2, sim_len - 1]],
            simulated_orbit[[3, sim_len - 1]]
        ),
    );
    new_system.array_parameters.insert(
        "vv0".to_string(),
        array!(
            simulated_orbit[[4, sim_len - 1]],
            simulated_orbit[[5, sim_len - 1]],
            simulated_orbit[[6, sim_len - 1]]
        ),
    );
    new_system.steps = 2;
    new_system
}

/// Simulates orbit onward from an impulse
pub fn simulate_orbit_from_impulse(
    prev_impulsed_orbit: Array2<f64>,
    system: &System,
    impulse_t: f64,
) -> Array2<f64> {
    //Add a patch point to the orbit if simulation start time is after the impulse
    let mut impulsed_system: System = system.clone();
    let start_t: f64 = system.start_time;
    let mut impulsed_orbit: Array2<f64> = prev_impulsed_orbit.clone();
    if (!impulse_t.is_nan()) && (impulse_t < start_t) {
        let impulse_point: Array2<f64> = array![
            [impulse_t],
            [system.array_parameters["rr0"][0]],
            [system.array_parameters["rr0"][1]],
            [system.array_parameters["rr0"][2]],
            [system.array_parameters["vv0"][0]],
            [system.array_parameters["vv0"][1]],
            [system.array_parameters["vv0"][2]]
        ];
        impulsed_orbit =
            concatenate(Axis(1), &[impulsed_orbit.view(), impulse_point.view()]).unwrap();
        let tmp_system =
            create_patch_system_orbit_to_impulse(system.clone(), start_t, &impulse_point);
        let patch_orbit: Array2<f64> = simulate_system(tmp_system);
        impulsed_system.array_parameters.insert(
            "rr0".to_string(),
            array!(
                patch_orbit[[1, 1]],
                patch_orbit[[2, 1]],
                patch_orbit[[3, 1]]
            ),
        );
        impulsed_system.array_parameters.insert(
            "vv0".to_string(),
            array!(
                patch_orbit[[4, 1]],
                patch_orbit[[5, 1]],
                patch_orbit[[6, 1]]
            ),
        );
    }

    let simulated_orbit: Array2<f64> = simulate_system(impulsed_system);
    //Add simulated patch to the final orbit
    if impulsed_orbit == array![[]] {
        impulsed_orbit = simulated_orbit.clone();
    } else {
        impulsed_orbit = concatenate(
            Axis(1),
            &[impulsed_orbit.view(), simulated_orbit.clone().view()],
        )
        .unwrap();
    }
    impulsed_orbit
}

/// Calculates the effect of a given impulse
///
/// Inputs:
///
/// system: a system struct TODO: add link within docs
///
/// orbit: the orbit unperturbed by the impulse TODO: which units
///
/// impulse: an array of changes in coordinates (TODO: which ones) due to the impulse at given times [time-array, coordinate-arrays, mass-change-array]
///
/// Output: An array containing the orbit modified by the impulse
pub fn calculate_impulse_effect(system: System, impulse: Array2<f64>) -> Array2<f64> {
    let system_step_size: f64 = (system.end_time - system.start_time) / (system.steps as f64 - 1.);
    let mut impulsed_orbit: Array2<f64> = array![[]];
    let mut prev_ti: f64 = f64::NAN;
    let mut mut_system: System = system.clone();
    impulse.row(0).into_iter().enumerate().for_each(|(i, ti)| {
        let steps: f64 = ((*ti - mut_system.start_time) / system_step_size).ceil();
        mut_system.end_time = mut_system.start_time + (steps - 1.) * system_step_size;
        mut_system.steps = steps as usize;
        impulsed_orbit = simulate_orbit_from_impulse(impulsed_orbit.clone(), &mut_system, prev_ti);

        //An intermediate system to simulate from last orbit time point to the impulse time point
        let tmp_system: System =
            create_patch_system_orbit_to_impulse(mut_system.clone(), *ti, &impulsed_orbit);
        let patch_orbit: Array2<f64> = simulate_system(tmp_system);

        //The new mu and initial conditions after applying the impulse at ti
        mut_system.f64_parameters.insert(
            "mu".to_string(),
            mut_system.f64_parameters["mu"] + impulse[[4, i]] / G,
        );
        mut_system.array_parameters.insert(
            "rr0".to_string(),
            array!(
                patch_orbit[[1, 1]],
                patch_orbit[[2, 1]],
                patch_orbit[[3, 1]]
            ),
        );
        mut_system.array_parameters.insert(
            "vv0".to_string(),
            array!(
                patch_orbit[[4, 1]],
                patch_orbit[[5, 1]],
                patch_orbit[[6, 1]]
            ) + array!(impulse[[1, i]], impulse[[2, i]], impulse[[3, i]]),
        );

        mut_system.start_time = mut_system.end_time + system_step_size;
        prev_ti = *ti;
    });
    let steps: f64 = ((system.end_time - mut_system.start_time) / system_step_size + 1.).round();
    mut_system.steps = steps as usize;
    mut_system.end_time = system.end_time;
    simulate_orbit_from_impulse(impulsed_orbit, &mut_system, prev_ti)
}

/// Calculates the effect of a given impulse without mass change on an orbit
///
/// **Inputs**:
///
/// orbit: the orbit unperturbed by the impulse TODO: which units
///
/// impulse: an array of changes in coordinates (TODO: which ones) due to the impulse at given times [time-array, coordinate-arrays]
///
/// **Output**: An array containing the orbit modified by the impulse
pub fn calculate_impulse_effect_without_mass_change(
    system: System,
    impulse: Array2<f64>,
) -> Array2<f64> {
    let mass_changes: Array2<f64> = Array2::zeros((1, impulse.row(0).len()));
    calculate_impulse_effect(
        system,
        concatenate(Axis(0), &[impulse.view(), mass_changes.view()]).unwrap(),
    )
}

//TODO: write tests
//TODO: calculate orbit shapes only
