use math::{cross_product, euclidean_norm};
use ndarray::Array1;
use std::f64::consts::PI;

pub mod math;

pub fn calculate_mu(m1: f64, m2: f64) -> f64 {
    let gravitational_constant: f64 = 6.67430 * (10_f64).powf(-11.);
    gravitational_constant * (m1 + m2)
}

pub fn calculate_eccentricity(position: Array1<f64>, velocity: Array1<f64>, mu: f64) -> f64 {
    let angular_momentum_per_unit_mass: Array1<f64> =
        cross_product(position.clone(), velocity.clone());
    let e: Array1<f64> = -cross_product(angular_momentum_per_unit_mass, velocity) / mu
        - position.clone() / euclidean_norm(position);
    euclidean_norm(e)
}

pub fn calculate_mean_anomaly(t: Array1<f64>, n: f64, tau: f64) -> Array1<f64> {
    n * (t - tau)
}

pub fn calculate_true_anomaly_from_series(
    t: Array1<f64>,
    eccentricity: f64,
    rotation_time: f64,
    periapsis_time: f64,
) -> Array1<f64> {
    let mean_anomaly = calculate_mean_anomaly(t, 2. * PI / rotation_time, periapsis_time);
    mean_anomaly.clone()
        + (2. * eccentricity - eccentricity.powf(3.) / 4.)
            * mean_anomaly.clone().mapv_into(|v| v.sin())
        + (5. * eccentricity.powf(2.) / 4.) * (2. * mean_anomaly.clone()).mapv_into(|v| v.sin())
        + (13. * eccentricity.powf(3.) / 12.) * (3. * mean_anomaly.clone()).mapv_into(|v| v.sin())
}

pub fn calculate_radius_from_true_anomaly(
    true_anomaly: Array1<f64>,
    eccentricity: f64,
    semimajor_axis: f64,
) -> Array1<f64> {
    semimajor_axis * (1. - eccentricity.powf(2.))
        / (1. + eccentricity * true_anomaly.mapv_into(|v| v.cos()))
}
