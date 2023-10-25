use ndarray::Array1;
use std::f64::consts::PI;

pub fn calculate_mean_anomaly(
    t: Array1<f64>,
    rotation_time: f64,
    periapsis_time: f64,
) -> Array1<f64> {
    2.0 * PI * (t - periapsis_time) / rotation_time
}

pub fn calculate_true_anomaly_from_series(
    t: Array1<f64>,
    eccentricity: f64,
    rotation_time: f64,
    periapsis_time: f64,
) -> Array1<f64> {
    let mean_anomaly = calculate_mean_anomaly(t, rotation_time, periapsis_time);
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
