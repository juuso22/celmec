use crate::system::System;
use ndarray::{array, concatenate, Array2, Axis};

/// Calculates the effect of a given impulse
///
/// **Inputs**:
///
/// **system**: a system struct TODO: add link within docs
///
/// **impulse_orig**: an array of changes in coordinates (TODO: which ones) due to the impulse at given times [time-array, coordinate-arrays, mass-change-array]
///
/// **Output**: An array containing the orbit modified by the impulse
pub fn calculate_impulse_effect(
    system: &(impl System + std::fmt::Debug),
    impulse_orig: Array2<f64>,
) -> Array2<f64> {
    let mut impulse: Array2<f64> = impulse_orig.clone();
    if (impulse.dim().0 != 4) | (impulse.dim().0 != 4) {
        panic!("Cannot interpret impulse of dimensions {:#?}. At least 4 rows are needed (t, v_x, v_y and v_z)", impulse.dim());
    } else if impulse.dim().0 == 4 {
        //If no mass ajection was given, add mass ejection equal to 0
        let mass_changes: Array2<f64> = Array2::zeros((1, impulse.row(0).len()));
        impulse = concatenate(Axis(0), &[impulse.view(), mass_changes.view()]).unwrap();
    }

    let mut impulsed_orbit: Array2<f64> = array![[]];
    let mut mut_system = system.clone();
    impulse.row(0).into_iter().enumerate().for_each(|(i, _ti)| {
        let steps: f64 =
            ((impulse[[0, i]] - mut_system.get_start_time()) / mut_system.get_step_size()).ceil();
        mut_system = mut_system
            .set_end_time(mut_system.get_start_time() + (steps - 1.) * mut_system.get_step_size());
        mut_system = mut_system.set_steps(steps as usize);

        let simulated_patch: Array2<f64> = mut_system.simulate();
        if i == 0 {
            impulsed_orbit = simulated_patch.clone();
        } else {
            impulsed_orbit = concatenate(
                Axis(1),
                &[impulsed_orbit.view(), simulated_patch.clone().view()],
            )
            .unwrap();
        }

        mut_system = mut_system.simulate_impulse(
            impulse.column(i).to_owned(),
            simulated_patch
                .column(simulated_patch.dim().1 - 1)
                .to_owned(),
        );
    });

    let steps: f64 =
        ((system.get_end_time() - mut_system.get_start_time()) / mut_system.get_step_size() + 1.)
            .round();
    mut_system = mut_system.set_steps(steps as usize);
    mut_system = mut_system.set_end_time(system.get_end_time());
    let simulated_patch: Array2<f64> = mut_system.simulate();
    concatenate(Axis(1), &[impulsed_orbit.view(), simulated_patch.view()]).unwrap()
}

//TODO: write tests
//TODO: calculate orbit shapes only
