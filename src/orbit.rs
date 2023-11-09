use math::{cross_product, euclidean_norm, solve_equation_iteratively};
use ndarray::Array1;
use std::collections::HashMap;
use std::f64::consts::PI;

pub mod math;

pub fn calculate_mu(m1: f64, m2: f64) -> f64 {
    let gravitational_constant: f64 = 6.67430 * (10_f64).powf(-11.);
    gravitational_constant * (m1 + m2)
}

pub fn calculate_e(position: Array1<f64>, velocity: Array1<f64>, mu: f64) -> Array1<f64> {
    let angular_momentum_per_unit_mass: Array1<f64> =
        cross_product(position.clone(), velocity.clone());
    -cross_product(angular_momentum_per_unit_mass, velocity) / mu
        - position.clone() / euclidean_norm(position)
}

pub fn calculate_eccentricity(position: Array1<f64>, velocity: Array1<f64>, mu: f64) -> f64 {
    euclidean_norm(calculate_e(position, velocity, mu))
}

pub fn calculate_h(position: Array1<f64>, velocity: Array1<f64>, mu: f64) -> f64 {
    0.5 * euclidean_norm(velocity).powf(2.) - mu / euclidean_norm(position)
}

pub fn calculate_a(mu: f64, h: f64, eccentricity: f64) -> f64 {
    let mut sign: f64 = -1.;
    if eccentricity > 1. {
        sign = 1.;
    }
    sign * mu / 2. / h
}

pub fn calculate_n(mu: f64, a: f64) -> f64 {
    mu.sqrt() * a.powf(-3. / 2.)
}

pub fn calculate_mean_anomaly(t: Array1<f64>, n: f64, periapsis_time: f64) -> Array1<f64> {
    n * (t - periapsis_time)
}

pub fn calculate_mean_anomaly_from_eccentric_anomaly(
    eccentric_anomaly: Array1<f64>,
    eccentricity: f64,
) -> Array1<f64> {
    if eccentricity > 1. {
        eccentric_anomaly.clone() - eccentricity * eccentric_anomaly.mapv_into(|v| v.sinh())
    } else if (eccentricity >= 0.) && (eccentricity < 1.) {
        eccentric_anomaly.clone() - eccentricity * eccentric_anomaly.mapv_into(|v| v.sin())
    } else if eccentricity == 1. {
        eccentric_anomaly.clone().mapv_into(|v| v.powf(3.)) / 6. + eccentric_anomaly / 2.
    } else {
        panic!("Eccentricity cannot be negative!")
    }
}

pub fn calculate_initial_true_anomaly_from_initial_conditions(
    position: Array1<f64>,
    e: Array1<f64>,
    eccentricity: f64,
) -> f64 {
    let r_len: f64 = euclidean_norm(position.clone());
    let initial_true_anomaly_cos: f64 = (position.clone() * e.clone()).sum() / r_len / eccentricity;
    // The following check implies sin of the true anomaly is negative, hence the angle is in the range (-PI, 0)
    if euclidean_norm(cross_product(position, e)) < 0. {
        return -initial_true_anomaly_cos.acos();
    }
    initial_true_anomaly_cos.acos()
}

pub fn calculate_eccentric_anomaly_from_true_anomaly(
    true_anomaly: Array1<f64>,
    eccentricity: f64,
) -> Array1<f64> {
    if eccentricity > 1. {
        let true_anomaly_cos: Array1<f64> = true_anomaly.mapv_into(|v| v.cos());
        let hyperbolic_anomaly_cosh: Array1<f64> =
            (true_anomaly_cos.clone() + eccentricity) / (1. + true_anomaly_cos * eccentricity);
        hyperbolic_anomaly_cosh.mapv_into(|v| v.acosh())
    } else if (eccentricity >= 0.) && (eccentricity < 1.) {
        let eccentric_anomaly_cos: Array1<(f64, f64)> = true_anomaly.mapv_into_any(|v| {
            (
                v.sin().signum(),
                (v.cos() + eccentricity) / (1. + eccentricity * v.cos()),
            )
        });
        eccentric_anomaly_cos.mapv_into_any(|v| v.0 * v.1.acos())
    } else if eccentricity == 1. {
        (true_anomaly / 2.).mapv_into(|v| v.tan())
    } else {
        panic!("Eccentricity cannot be negative!")
    }
}

pub fn calculate_periapsis_time(t: f64, mean_anomaly: f64, n: f64) -> f64 {
    t - mean_anomaly / n
}

fn kepler_eq_iterative_step(
    eccentric_anomaly: Array1<f64>,
    time: Array1<f64>,
    parameters: HashMap<&str, f64>,
) -> Array1<f64> {
    parameters["eccentricity"] * eccentric_anomaly.mapv_into(|v| v.sin())
        + parameters["n"] * (time - parameters["periapsis_time"])
}

fn hyperbolic_kepler_eq_iterative_step(
    hyperbolic_anomaly: Array1<f64>,
    time: Array1<f64>,
    parameters: HashMap<&str, f64>,
) -> Array1<f64> {
    let mean_anomaly: Array1<f64> = parameters["n"] * (time - parameters["periapsis_time"]);
    let f: Array1<f64> = parameters["eccentricity"]
        * hyperbolic_anomaly.clone().mapv_into(|v| v.sinh())
        - hyperbolic_anomaly.clone()
        - mean_anomaly;
    let df_dh: Array1<f64> =
        parameters["eccentricity"] * hyperbolic_anomaly.clone().mapv_into(|v| v.cosh()) - 1.;
    hyperbolic_anomaly - f / df_dh
    //    parameters["eccentricity"] * eccentric_anomaly.mapv_into(|v| v.sinh())
    //        + parameters["n"] * (time - parameters["periapsis_time"])
}

fn barker_eq_iterative_step(
    eccentric_anomaly: Array1<f64>,
    time: Array1<f64>,
    parameters: HashMap<&str, f64>,
) -> Array1<f64> {
    2. * parameters["n"] * (time - parameters["periapsis_time"])
        - eccentric_anomaly.mapv_into(|v| v.powf(3.)) / 3.
}

pub fn calculate_eccentric_anomaly_from_kepler_equation(
    t: Array1<f64>,
    initial_value: Array1<f64>,
    tolerance: f64,
    max_iterations: usize,
    n: f64,
    eccentricity: f64,
    periapsis_time: f64,
) -> Array1<f64> {
    let parameters: HashMap<&str, f64> = HashMap::from([
        ("n", n),
        ("eccentricity", eccentricity),
        ("periapsis_time", periapsis_time),
    ]);
    if eccentricity > 1. {
        solve_equation_iteratively(
            &hyperbolic_kepler_eq_iterative_step,
            initial_value,
            t,
            parameters,
            tolerance,
            max_iterations,
        )
    } else if eccentricity == 1. {
        solve_equation_iteratively(
            &barker_eq_iterative_step,
            initial_value,
            t,
            parameters,
            tolerance,
            max_iterations,
        )
    } else if (eccentricity < 1.) && (eccentricity >= 0.) {
        //        let initial_value: Array1<f64> = t.clone().mapv_into(|v| {
        //            eccentricity * (n * v - periapsis_time).sin() / (1. - eccentricity * (n * (v - periapsis_time)).cos())
        //        });
        solve_equation_iteratively(
            &kepler_eq_iterative_step,
            initial_value,
            t,
            parameters,
            tolerance,
            max_iterations,
        )
    } else {
        panic!("Eccentricity cannot be negative.");
    }
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

pub fn calculate_true_anomaly_from_eccentric_anomaly(
    eccentric_anomaly: Array1<f64>,
    eccentricity: f64,
) -> Array1<f64> {
    if eccentricity > 1. {
        (((eccentricity + 1.) / (eccentricity - 1.)).sqrt()
            * (eccentric_anomaly / 2.).mapv_into(|v| v.tanh()))
        .mapv_into(|v| v.atan())
            * 2.
    } else if (eccentricity < 1.) && (eccentricity >= 0.) {
        let true_anomaly_cos: Array1<(f64, f64)> = eccentric_anomaly.mapv_into_any(|v| {
            (
                v.sin().signum(),
                (v.cos() - eccentricity) / (1. - eccentricity * v.cos()),
            )
        });
        true_anomaly_cos.mapv_into_any(|v| v.0 * v.1.acos())
    } else if eccentricity == 1. {
        eccentric_anomaly.mapv_into(|v| v.atan()) / 2.
    } else {
        panic!("Eccentricity cannot be negative!")
    }
}

pub fn calculate_radius_from_true_anomaly(
    true_anomaly: Array1<f64>,
    eccentricity: f64,
    semimajor_axis: f64,
) -> Array1<f64> {
    semimajor_axis * (1. - eccentricity.powf(2.)).abs()
        / (1. + eccentricity * true_anomaly.mapv_into(|v| v.cos()))
}
