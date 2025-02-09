use crate::orbital_elements::{
    calculate_keplerian_elements_from_initial_rr_and_vv_and_mu, KeplerianElements,
};
use crate::system::System;
use crate::transformations::cartesian_coordinates_from_f_r_and_keplerian_elements;
use crate::two_body;
use ndarray::{concatenate, Array1, Array2, Axis};

/// Does a simulation based on a system
///
/// **Inputs**:
///
/// system: a system struct with info about the system to be simulated
///
/// **Output**: An array of arrays containing the simulated orbit. The inner arrays are, in this order: t, x, y, z, v_x, v_y, v_z
pub fn simulate_system(system: System) -> Array2<f64> {
    let t: Array1<f64> = Array1::linspace(system.start_time, system.end_time, system.steps);
    if system.model == "2-body".to_string() {
        //TODO: have each model file implement a trait (to be added to this file) that does the simulation ie. what comes below for 2-body
        if system.perturbations.is_empty() {
            let keplerian_elements: KeplerianElements =
                calculate_keplerian_elements_from_initial_rr_and_vv_and_mu(
                    system.array_parameters["rr0"].clone(),
                    system.array_parameters["vv0"].clone(),
                    system.f64_parameters["mu"],
                    system.start_time,
                );
            let eccentric_anomaly: Array1<f64> =
                two_body::calculate_eccentric_anomaly_from_initial_rr_and_vv(
                    system.array_parameters["rr0"].clone(),
                    system.array_parameters["vv0"].clone(),
                    system.f64_parameters["mu"],
                    system.start_time,
                    system.end_time,
                    system.steps,
                );
            let f: Array1<f64> = two_body::calculate_f_from_eccentric_anomaly(
                eccentric_anomaly.clone(),
                keplerian_elements.e,
            );
            let r: Array1<f64> =
                two_body::calculate_r_from_f(f.clone(), keplerian_elements.e, keplerian_elements.a);
            let v: Array1<f64> = two_body::calculate_v_from_eccentric_anomaly(
                eccentric_anomaly,
                system.f64_parameters["mu"],
                keplerian_elements.e,
                keplerian_elements.a,
            );
            let rr: Array2<f64> =
                cartesian_coordinates_from_f_r_and_keplerian_elements(f, r, keplerian_elements);
            let vv: Array2<f64> = two_body::calculate_vv_from_v_rr_and_initial_vv(
                v,
                rr.clone(),
                system.array_parameters["vv0"].clone(),
            );
            let mut t_2d: Array2<f64> = Array2::zeros((1, t.len()));
            t_2d.row_mut(0).assign(&t);
            concatenate(Axis(0), &[t_2d.view(), rr.view(), vv.view()]).unwrap()
        } else {
            panic!(
                "Orbit simulation for model {} with perturbations has not been implemented",
                system.model
            );
        }
    } else {
        panic!(
            "Orbit simulation for model {} has not been implemented",
            system.model
        );
    }
}
