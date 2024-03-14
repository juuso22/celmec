use math::{cross_product, euclidean_norm, solve_equation_iteratively};
use ndarray::Array1;
use std::collections::HashMap;
use std::f64::consts::PI;

pub mod math;

pub fn calculate_mu(m1: f64, m2: f64) -> f64 {
    let gravitational_constant: f64 = 6.67430 * (10_f64).powf(-11.);
    gravitational_constant * (m1 + m2)
}

pub fn calculate_ee(rr: Array1<f64>, vv: Array1<f64>, mu: f64) -> Array1<f64> {
    let angular_momentum_per_unit_mass: Array1<f64> =
        cross_product(rr.clone(), vv.clone());
    -cross_product(angular_momentum_per_unit_mass, vv) / mu
        - rr.clone() / euclidean_norm(rr)
}

pub fn calculate_e(rr: Array1<f64>, vv: Array1<f64>, mu: f64) -> f64 {
    euclidean_norm(calculate_ee(rr, vv, mu))
}

pub fn calculate_h(rr: Array1<f64>, vv: Array1<f64>, mu: f64) -> f64 {
    0.5 * euclidean_norm(vv).powf(2.) - mu / euclidean_norm(rr)
}

pub fn calculate_a(mu: f64, h: f64, e: f64) -> f64 {
    let mut sign: f64 = -1.;
    if e > 1. {
        sign = 1.;
    }
    sign * mu / 2. / h
}

pub fn calculate_n(mu: f64, a: f64) -> f64 {
    mu.sqrt() * a.powf(-3. / 2.)
}

pub fn calculate_mean_anomaly(t: Array1<f64>, n: f64, tau: f64) -> Array1<f64> {
    n * (t - tau)
}

pub fn calculate_mean_anomaly_from_eccentric_anomaly(
    eccentric_anomaly: Array1<f64>,
    e: f64,
) -> Array1<f64> {
    if e > 1. {
        eccentric_anomaly.clone() - e * eccentric_anomaly.mapv_into(|v| v.sinh())
    } else if (e >= 0.) && (e < 1.) {
        eccentric_anomaly.clone() - e * eccentric_anomaly.mapv_into(|v| v.sin())
    } else if e == 1. {
        eccentric_anomaly.clone().mapv_into(|v| v.powf(3.)) / 6. + eccentric_anomaly / 2.
    } else {
        panic!("Eccentricity cannot be negative!")
    }
}

pub fn calculate_initial_f_from_initial_conditions(
    rr: Array1<f64>,
    ee: Array1<f64>,
    e: f64,
) -> f64 {
    let r: f64 = euclidean_norm(rr.clone());
    let initial_f_cos: f64 =
        (rr.clone() * ee.clone()).sum() / r / e;
    // The following check implies sin of the true anomaly is negative, hence the angle is in the range (-PI, 0)
    if euclidean_norm(cross_product(rr, ee)) < 0. {
        return -initial_f_cos.acos();
    }
    initial_f_cos.acos()
}

pub fn calculate_eccentric_anomaly_from_f(
    f: Array1<f64>,
    e: f64,
) -> Array1<f64> {
    if e > 1. {
        let f_cos: Array1<f64> = f.mapv_into(|v| v.cos());
        let hyperbolic_anomaly_cosh: Array1<f64> =
            (f_cos.clone() + e) / (1. + f_cos * e);
        hyperbolic_anomaly_cosh.mapv_into(|v| v.acosh())
    } else if (e >= 0.) && (e < 1.) {
        let eccentric_anomaly_cos: Array1<(f64, f64)> = f.mapv_into_any(|v| {
            (
                v.sin().signum(),
                (v.cos() + e) / (1. + e * v.cos()),
            )
        });
        eccentric_anomaly_cos.mapv_into_any(|v| v.0 * v.1.acos())
    } else if e == 1. {
        (f / 2.).mapv_into(|v| v.tan())
    } else {
        panic!("Eccentricity cannot be negative!")
    }
}

pub fn calculate_tau(t: f64, mean_anomaly: f64, n: f64) -> f64 {
    t - mean_anomaly / n
}

fn kepler_eq_iterative_step(
    eccentric_anomaly: Array1<f64>,
    time: Array1<f64>,
    parameters: HashMap<&str, f64>,
) -> Array1<f64> {
    parameters["e"] * eccentric_anomaly.mapv_into(|v| v.sin())
        + parameters["n"] * (time - parameters["tau"])
}

fn hyperbolic_kepler_eq_iterative_step(
    hyperbolic_anomaly: Array1<f64>,
    time: Array1<f64>,
    parameters: HashMap<&str, f64>,
) -> Array1<f64> {
    let mean_anomaly: Array1<f64> = parameters["n"] * (time - parameters["tau"]);
    let f: Array1<f64> = parameters["e"]
        * hyperbolic_anomaly.clone().mapv_into(|v| v.sinh())
        - hyperbolic_anomaly.clone()
        - mean_anomaly;
    let df_dh: Array1<f64> =
        parameters["e"] * hyperbolic_anomaly.clone().mapv_into(|v| v.cosh()) - 1.;
    hyperbolic_anomaly - f / df_dh
    //    parameters["e"] * eccentric_anomaly.mapv_into(|v| v.sinh())
    //        + parameters["n"] * (time - parameters["tau"])
}

fn barker_eq_iterative_step(
    eccentric_anomaly: Array1<f64>,
    time: Array1<f64>,
    parameters: HashMap<&str, f64>,
) -> Array1<f64> {
    2. * parameters["n"] * (time - parameters["tau"])
        - eccentric_anomaly.mapv_into(|v| v.powf(3.)) / 3.
}

pub fn calculate_eccentric_anomaly_iteratively(
    t: Array1<f64>,
    initial_value: Array1<f64>,
    tolerance: f64,
    max_iterations: usize,
    n: f64,
    e: f64,
    tau: f64,
) -> Array1<f64> {
    let parameters: HashMap<&str, f64> = HashMap::from([
        ("n", n),
        ("e", e),
        ("tau", tau),
    ]);
    if e > 1. {
        solve_equation_iteratively(
            &hyperbolic_kepler_eq_iterative_step,
            initial_value,
            t,
            parameters,
            tolerance,
            max_iterations,
        )
    } else if e == 1. {
        solve_equation_iteratively(
            &barker_eq_iterative_step,
            initial_value,
            t,
            parameters,
            tolerance,
            max_iterations,
        )
    } else if (e < 1.) && (e >= 0.) {
        //        let initial_value: Array1<f64> = t.clone().mapv_into(|v| {
        //            e * (n * v - tau).sin() / (1. - e * (n * (v - tau)).cos())
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

pub fn calculate_f_from_series(
    t: Array1<f64>,
    e: f64,
    rotation_time: f64,
    tau: f64,
) -> Array1<f64> {
    let mean_anomaly = calculate_mean_anomaly(t, 2. * PI / rotation_time, tau);
    mean_anomaly.clone()
        + (2. * e - e.powf(3.) / 4.)
            * mean_anomaly.clone().mapv_into(|v| v.sin())
        + (5. * e.powf(2.) / 4.) * (2. * mean_anomaly.clone()).mapv_into(|v| v.sin())
        + (13. * e.powf(3.) / 12.) * (3. * mean_anomaly.clone()).mapv_into(|v| v.sin())
}

pub fn calculate_f_from_eccentric_anomaly(
    eccentric_anomaly: Array1<f64>,
    e: f64,
) -> Array1<f64> {
    if e > 1. {
        (((e + 1.) / (e - 1.)).sqrt()
            * (eccentric_anomaly / 2.).mapv_into(|v| v.tanh()))
        .mapv_into(|v| v.atan())
            * 2.
    } else if (e < 1.) && (e >= 0.) {
        let f_cos: Array1<(f64, f64)> = eccentric_anomaly.mapv_into_any(|v| {
            (
                v.sin().signum(),
                (v.cos() - e) / (1. - e * v.cos()),
            )
        });
        f_cos.mapv_into_any(|v| v.0 * v.1.acos())
    } else if e == 1. {
        eccentric_anomaly.mapv_into(|v| v.atan()) / 2.
    } else {
        panic!("Eccentricity cannot be negative!")
    }
}

pub fn calculate_radius_from_f(
    f: Array1<f64>,
    e: f64,
    a: f64,
) -> Array1<f64> {
    a * (1. - e.powf(2.)).abs()
        / (1. + e * f.mapv_into(|v| v.cos()))
}
