use crate::constants::G;
use crate::math::{cross_product, euclidean_norm, solve_equation_iteratively};
use crate::orbital_elements::{
    calculate_a, calculate_e, calculate_ee,
    calculate_keplerian_elements_from_initial_rr_and_vv_and_mu, KeplerianElements,
};
use crate::transformations::cartesian_coordinates_from_f_r_and_keplerian_elements;
use ndarray::{array, Array1, Array2};
use std::collections::HashMap;
use std::f64::consts::PI;

/// Calculates μ (mu) from two masses m1 and m2.
///
/// μ = (m1 + m2) * [`G`],
///
/// where
///
/// G = gravitational constant
///
/// μ can be used to calculate the gravitational force between the 2 masses if the distance r between them is known:
///
/// ```
/// use celmec::orbit::calculate_mu;
///
/// let m1: f64 = 10.;
/// let m2: f64 = 20.;
/// let r: f64 = 100.;
/// let force = calculate_mu(m1, m2) / r.powf(2.);
/// ```
pub fn calculate_mu(m1: f64, m2: f64) -> f64 {
    G * (m1 + m2)
}

/// Calculates the Lagrangian h of a 2-body system.
///
/// **Inputs**:
///
/// rr = position of the 2 bodies with respect to each other
///
/// vv = velocity of the 2 bodies with respect to each other
///
/// μ = see [μ](`calculate_mu`).
///
/// **Output**: Lagrangian h
///
/// For a Newtonian 2-body system the Lagrangian is the difference between the kinetic and the potential energy of the system:
///
/// h = 0.5 * |**v**|<sup>2</sup> - μ / |**r**|,
///
/// **r** = position of the 2 bodies relative to each other at some point in time,
///
/// **v** = velocity of the 2 bodies relative to each other at the same point of time as **r**
///
/// μ = see [μ](`calculate_mu`)
///
/// Inputs are a position **r** (`rr`) and the velocity **v** (`vv`) of the 2 bodies with respect to each other at any one point of time as well as their [μ](`calculate_mu`).
pub fn calculate_h(rr: Array1<f64>, vv: Array1<f64>, mu: f64) -> f64 {
    0.5 * euclidean_norm(vv).powf(2.) - mu / euclidean_norm(rr)
}

/// Calculates focal parameter for parabola from a known distance and eccentric anomaly at some point of time for a 2-body system
///
/// **Inputs**:
///
/// eccentric_anomaly: eccentric anomaly
///
/// r: distance
///
/// **Output**: focal parameter of a parabola
///
/// The relevant equation is:
///
/// p = 2*r / (E<sup>2</sup> + 1),
///
/// where
///
/// p = the focal parameter of a conic section
///
/// r = distance
///
/// E = eccentric anomaly
pub fn calculate_focal_parameter_for_parabola_from_orbital_position(
    eccentric_anomaly: f64,
    r: f64,
) -> f64 {
    2. * r / (eccentric_anomaly.powf(2.0) + 1.)
}

/// Calculates the average angular velocity in a 2-body system
///
/// **Inputs**
///
/// mu: see [μ](`calculate_mu`)
///
/// a: semi-major axis [a](`calculate_a`) or [focal parameter p](`calculate_focal_parameter_for_parabola_from_orbital_position`) (in case of a parabolic orbit)
///
/// **Output**: Average angular velocity
///
/// n = μ<sup>1/2</sup> * a<sup>-3/2</sup>
///
/// where
///
/// n = average angular velocity
///
/// μ = see [μ](`calculate_mu`)
///
/// a = Semi-major axis or [focal parameter](`calculate_focal_parameter_for_parabola_from_orbital_position`) (in case of a parabolic orbit)
///
/// For elliptic orbits this is the same as calculating 2π / P, where P is the period of the orbit.
pub fn calculate_n(mu: f64, a: f64) -> f64 {
    mu.sqrt() * a.powf(-3. / 2.)
}

/// Calculates the mean anomaly M in a 2-body system
///
/// M = n * (t - τ)
///
/// where
///
/// n = average angular velocity
///
/// t = time
///
/// τ = perihelion time
///
/// Inputs are time t as well as average angular velocity [n](`calculate_n`) and perihelion time τ of the system.
///
/// If mean anomaly at some point of time is known, the equation above can be used to [calculate τ](`calculate_tau`).
pub fn calculate_mean_anomaly(t: Array1<f64>, n: f64, tau: f64) -> Array1<f64> {
    n * (t - tau)
}

/// calculates mean anomaly M from eccentric anomaly E for a 2-body system.
///
/// **Inputs**:
///
/// eccentric_anomaly: an array of eccentric anomalies
///
/// e: eccentricity
///
/// **Output**: an array of true anomalies
///
/// A different equations is used depending on the value of eccentricity:
///
/// 0 <= e < 1: Kepler's equation M = E - e * sin(E)
///
/// e = 1: Barker's equation M = E<sup>3</sup> / 6 + E / 2
///
/// e > 1: 'hyperbolic Kepler's equation' M = E - e * sinh(E)
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

/// Calculates true anomaly f from the initial relative distance, the [`ee`](`calculate_ee`) vector and [eccentricity](`calculate_e`) of a 2-body system.
///
/// **Inputs**:
///
/// rr: distance between the two bodies at some point of time
///
/// ee: see [`ee`](`calculate_ee`)
///
/// e: eccentricity
///
/// **Output**: true anomaly at the position `rr`
///
/// f calculated from the formula:
///
/// f = cos(**r** &middot **e**) / (r * e),
///
/// where
///
/// **r** = `rr`
///
/// **e** = `ee`
///
/// r = |**r**|
///
/// e = eccentricity
///
/// Lenght of the cross product of **r** and **e** is used to determine sin(f) which in term is used to determine whether f is in range [-π, 0] or (0, π].
pub fn calculate_initial_f_from_initial_conditions(
    rr: Array1<f64>,
    ee: Array1<f64>,
    e: f64,
) -> f64 {
    let r: f64 = euclidean_norm(rr.clone());
    let initial_f_cos: f64 = (rr.clone() * ee.clone()).sum() / r / e;
    // The following check implies sin of the true anomaly is negative, hence the angle is in the range (-PI, 0)
    if euclidean_norm(cross_product(rr, ee)) < 0. {
        return -initial_f_cos.acos();
    }
    initial_f_cos.acos()
}

/// Calculates eccentric anomaly E from from true anomaly f and eccenctricity e for a 2-body system
///
/// **Inputs**:
///
/// f: array of true anomalies
///
/// e: [eccentricity](`calculate_e`)
///
/// **Output**: an array of eccentric anomalies
///
/// The equation to solve depends on which shape (elliptic, hyperbolic, parabolic) of the orbit ie. the value of e:
///
/// 0 <= e < 1 (elliptic orbit): E = acos((cos(f) + e) / (1 + e * cos(f)))
///
/// e = 1 (parabolic orbit): E = tan(f / 2)
///
/// e > 1 (hyperbolic orbit): E = acosh((cos(f) + e) / (1 + cos(f) * e))
///
/// For the elliptic and hyperbolic cases, sin(f) is used whether E is in range [-π, 0] or (0, π].
pub fn calculate_eccentric_anomaly_from_f(f: Array1<f64>, e: f64) -> Array1<f64> {
    if e > 1. {
        let f_cos: Array1<f64> = f.mapv_into(|v| v.cos());
        let eccentric_anomaly_cosh: Array1<f64> = (f_cos.clone() + e) / (1. + f_cos * e);
        eccentric_anomaly_cosh.mapv_into(|v| v.acosh())
    } else if (e >= 0.) && (e < 1.) {
        let eccentric_anomaly_cos: Array1<(f64, f64)> =
            f.mapv_into_any(|v| (v.sin().signum(), (v.cos() + e) / (1. + e * v.cos())));
        eccentric_anomaly_cos.mapv_into_any(|v| v.0 * v.1.acos())
    } else if e == 1. {
        (f / 2.).mapv_into(|v| v.tan())
    } else {
        panic!("Eccentricity cannot be negative!")
    }
}

/// Calculates perihelion time of a 2-body system
///
/// τ = t - M / n
///
/// where
///
/// t = time
///
/// n = average angular velocity
///
/// M = mean anomaly
///
/// Inputs are time t and the mean anomaly at that time as well as average angular velocity [n](`calculate_n`) of the system.
///
/// If τ is known, the equation above can be used to [calculate n](`calculate_n`).
pub fn calculate_tau(t: f64, mean_anomaly: f64, n: f64) -> f64 {
    t - mean_anomaly / n
}

/// Calculates on iterative step when solving eccentric anomaly from the the Kepler equation for a 2-body system.
///
/// The iterative step is:
///
/// E<sub>i+1</sub> = e * sin(E<sub>i</sub>) + n * (t - τ)
///
/// for the Kepler equation:
///
/// E - e*sin(E) = n * (t - τ),
///
/// where
///
/// E = eccentric anomaly
///
/// e = eccentricity
///
/// n = average angular velocity
///
/// t = time
///
/// τ = perihelion time
///
/// which cannot be solved in closed form for E.
///
/// Inputs are the eccentric anomalies at step i for some time range, the said time range and a hashmap of parameters that should contain n, e and τ.
///
/// The contents of the paremeters hashmap are currently not checked in any way so it is the responsibility of the user to make sure the right content is included.
///
/// Note that n * (t - τ) = M, where M is the mean anomaly.
fn kepler_eq_iterative_step(
    eccentric_anomaly: Array1<f64>,
    time: Array1<f64>,
    parameters: HashMap<&str, f64>,
) -> Array1<f64> {
    parameters["e"] * eccentric_anomaly.mapv_into(|v| v.sin())
        + parameters["n"] * (time - parameters["tau"])
}

/// Calculates on iterative step when solving eccentric anomaly E from the 'hyperbolic Kepler equation' for a 2-body system.
///
/// **Inputs**:
///
/// eccentric_anomaly: E at step i for some time range
///
/// time: the said time range. Note that the length of time and eccentric_anomaly should be the same
///
/// parameters: a hashmap of parameters that should contain n, e and τ (see their meaning below). The contents of the paremeters hashmap are currently not checked in any way so it is the responsibility of the user to make sure the right content is included.
///
/// **Output**: An array of eccentric anomalies.
///
/// The iterative step is:
///
/// E<sub>i+1</sub> = E<sub>i</sub> - (e * sinh(E<sub>i</sub>) - E<sub>i</sub> - n * (t - τ)) / (e * cosh(E<sub>i</sub>) - 1)
///
/// which is derived using Newton-Raphson method from the 'hyperbolic Kepler equation':
///
/// E - e*sinh(E) = n * (t - τ),
///
/// where
///
/// E = eccentric anomaly
///
/// e = eccentricity
///
/// n = see [n](`calculate_n`)
///
/// t = time
///
/// τ = perihelion time
///
/// which cannot be solved in closed form for E.
///
/// Note that n * (t - τ) = M, where M is the mean anomaly.
fn hyperbolic_kepler_eq_iterative_step(
    eccentric_anomaly: Array1<f64>,
    time: Array1<f64>,
    parameters: HashMap<&str, f64>,
) -> Array1<f64> {
    let mean_anomaly: Array1<f64> =
        calculate_mean_anomaly(time, parameters["n"], parameters["tau"]);
    let f: Array1<f64> = parameters["e"] * eccentric_anomaly.clone().mapv_into(|v| v.sinh())
        - eccentric_anomaly.clone()
        - mean_anomaly;
    let df_dh: Array1<f64> =
        parameters["e"] * eccentric_anomaly.clone().mapv_into(|v| v.cosh()) - 1.;
    eccentric_anomaly - f / df_dh
}

/// Calculates on iterative step when solving eccentric anomaly E from the Barker equation for a 2-body system.
///
/// **Inputs**:
///
/// eccentric_anomaly: E at step i for some time range
///
/// time: the said time range. Note that the length of time and eccentric_anomaly should be the same
///
/// parameters: a hashmap of parameters that should contain n and τ (see their meaning below). The contents of the paremeters hashmap are currently not checked in any way so it is the responsibility of the user to make sure the right content is included.
///
/// **Output**: An array of eccentric anomalies.
///
/// The iterative step is:
///
/// E<sub>i+1</sub> = E <sub>i</sub> - (E<sub>i</sub><sup>3</sup>/6 + E<sub>i</sub>/2 - n * (t - τ)) / (E<sub>i</sub><sup>2</sup>/3 + 1/2)
///
/// which is derived using Newton-Raphson method from the Barker equation:
///
/// E<sup>3</sup> / 6 + E / 2 = n * (t - τ),
///
/// where
///
/// E = eccentric anomaly
///
/// e = eccentricity
///
/// n = see [n](`calculate_n`)
///
/// t = time
///
/// τ = perihelion time
///
/// Note that n * (t - τ) = M, where M is the mean anomaly.
fn barker_eq_iterative_step(
    eccentric_anomaly: Array1<f64>,
    time: Array1<f64>,
    parameters: HashMap<&str, f64>,
) -> Array1<f64> {
    let f: Array1<f64> = eccentric_anomaly.clone().mapv_into(|v| v.powf(3.)) / 6.
        + eccentric_anomaly.clone() / 2.
        - parameters["n"] * (time - parameters["tau"]);
    let df_de: Array1<f64> = eccentric_anomaly.clone().mapv_into(|v| v.powf(2.)) / 3. + 0.5;
    eccentric_anomaly - f / df_de
}

/// Calculates the eccentric anomaly iteratively from some initial conditions.
///
/// **Inputs**:
///
/// t: time range for which we calculate
///
/// initial_value: initial guesses of the eccentric anomalies for the time range t
///
/// tolerance: stop iteration when the difference between the old and new value of the eccentric anomaly at a given time point is less than this parameter.
///
/// max_iterations: maximum number of itarations if tolerance is not reached first.
///
/// n: average angular velocity
///
/// e: eccentricity
///
/// tau: perihelion time
///
/// **Output**: An array of eccentric anomalies
///
/// The equation to be solved depends on the value of the eccentricity which is one of the inputs:
///
/// 0 <= e < 1: Kepler equation
///
/// e = 1: Barker equation
///
/// e > 1: 'hyperbolic Kepler equation'
pub fn calculate_eccentric_anomaly_iteratively(
    t: Array1<f64>,
    initial_value: Array1<f64>,
    tolerance: f64,
    max_iterations: usize,
    n: f64,
    e: f64,
    tau: f64,
) -> Array1<f64> {
    let parameters: HashMap<&str, f64> = HashMap::from([("n", n), ("e", e), ("tau", tau)]);
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

/// Calculates true anomaly f from Fourier series for a 2-body system.
///
/// **Warning**: this function is currently suitable only for elliptic orbits ie. 0 <= eccentricity e < 1.
///
/// **Inputs**:
///
/// t: a time range for which the true anomalies are solved
///
/// e: [eccentricity](`calculate_e`)
///
/// rotation_time: time of one rotation of one body around the other
///
/// tau: [perihelion time](`calculate_tau`) of the system
///
/// **Output**: An array of true anomalies.
pub fn calculate_f_from_series(
    t: Array1<f64>,
    e: f64,
    rotation_time: f64,
    tau: f64,
) -> Array1<f64> {
    let mean_anomaly = calculate_mean_anomaly(t, 2. * PI / rotation_time, tau);
    mean_anomaly.clone()
        + (2. * e - e.powf(3.) / 4.) * mean_anomaly.clone().mapv_into(|v| v.sin())
        + (5. * e.powf(2.) / 4.) * (2. * mean_anomaly.clone()).mapv_into(|v| v.sin())
        + (13. * e.powf(3.) / 12.) * (3. * mean_anomaly.clone()).mapv_into(|v| v.sin())
}

/// Calculates true anomaly f from eccentric anomaly E for a 2-body system.
///
/// **Inputs**:
///
/// eccentric_anomaly: an array of eccentric anomalies
///
/// e: [eccentricity](`calculate_e`)
///
/// rotation_time: time of one rotation of one body around the other
///
/// **Output**:
///
/// An array of true anomalies.
///
/// Diffrent equations are used for different orbit shapes ie. different eccentricites.
///
/// 0 <= e < 1: f = acos((cos(E) - e) / (1 - e * cos(E)))
///
/// e = 1:  f = atan(E) / 2
///
/// e > 1: f = atan(((e + 1) / (e - 1))<sup>1/2</sup> * tanh(E / 2))
pub fn calculate_f_from_eccentric_anomaly(eccentric_anomaly: Array1<f64>, e: f64) -> Array1<f64> {
    if e > 1. {
        (((e + 1.) / (e - 1.)).sqrt() * (eccentric_anomaly / 2.).mapv_into(|v| v.tanh()))
            .mapv_into(|v| v.atan())
            * 2.
    } else if (e < 1.) && (e >= 0.) {
        let f_cos: Array1<(f64, f64)> = eccentric_anomaly
            .mapv_into_any(|v| (v.sin().signum(), (v.cos() - e) / (1. - e * v.cos())));
        //sin(E) is used to determine whether f is in [-PI, 0] or [0, PI]
        f_cos.mapv_into_any(|v| v.0 * v.1.acos())
    } else if e == 1. {
        let f_cos: Array1<(f64, f64)> =
            eccentric_anomaly.mapv_into_any(|v| (v.sin().signum(), 2. / (v.powf(2.) + 1.) - 1.));
        //sin(E) is used to determine whether f is in [-PI, 0] or [0, PI]
        f_cos.mapv_into_any(|v| v.0 * v.1.acos())
    } else {
        panic!("Eccentricity cannot be negative!")
    }
}

/// Calculates the distance r between 2 bodies from true anomaly f, eccentricity e and semi-major axis a.
///
/// **Inputs**:
///
/// f: An array of true anomalies
///
/// e: eccentricity
///
/// a: semi-major axis
///
/// **Output**: An array of radii between to bodies attracted by gravitation
///
/// Using a reference frame with one of the bodies fixed as origin, this can also be undrstood as the radius of the orbit of the other body with respect to the first.
///
/// r = a * |(1-e<sup>2</sup>>)| / (1 + e * cos(f))
///
/// where
///
/// a = semi-major axis
///
/// e = eccentricity
///
/// f = true anomaly
///
/// Inputs are an array of true anomalies as well the semi-major axis [a](`calculate_a`) and the eccentricity [e](`calculate_e`) of the system.
pub fn calculate_r_from_f(f: Array1<f64>, e: f64, a: f64) -> Array1<f64> {
    if (e != 1.) && (e >= 0.) {
        a * (1. - e.powf(2.)).abs() / (1. + e * f.mapv_into(|v| v.cos()))
    } else if e == 1. {
        0.5 * a / (1. + f.mapv_into(|v| v.cos()))
    } else {
        panic!("Eccentricity cannot be negative!")
    }
}

/// Calculates the true anomaly f for a 2 body problem from initial conditions for a given total time split into a given number of steps.
///
/// One of the bodies lies at the origin and all the inputs are given with respect to this body, referred to as "the central body". The other body is referred to as "the rotating body"
///
/// **Inputs**:
///
/// rr: The initial position of the rotating body with respect to the central body.
///
/// vv: The initial velocity of the rotating body with respect to the central body.
///
/// mu: The gravitational parameter of the system. See [μ](`calculate_mu`).
///
/// total_time: The time to be simulated.
///
/// steps: The number of intervals the original interval [0, total_time] is split into.
///
/// **Output**: An array of true anomalies.
pub fn calculate_f_from_initial_rr_and_vv(
    rr: Array1<f64>,
    vv: Array1<f64>,
    mu: f64,
    start_time: f64,
    end_time: f64,
    steps: usize,
) -> Array1<f64> {
    let ee: Array1<f64> = calculate_ee(rr.clone(), vv.clone(), mu);
    let e: f64 = calculate_e(rr.clone(), vv.clone(), mu);
    let initial_f: f64 = calculate_initial_f_from_initial_conditions(rr.clone(), ee, e);
    let initial_eccentric_anomaly: Array1<f64> =
        calculate_eccentric_anomaly_from_f(array![initial_f], e);
    let a: f64 = calculate_a_from_initial_rr_and_vv(rr, vv, mu);
    let n: f64 = calculate_n(mu, a);
    let tau: f64 = calculate_tau(
        0.,
        calculate_mean_anomaly_from_eccentric_anomaly(initial_eccentric_anomaly, e)[0],
        n,
    );

    let t: Array1<f64> = Array1::linspace(start_time, end_time, steps);
    let eccentric_anomaly: Array1<f64> = calculate_eccentric_anomaly_iteratively(
        t.clone(),
        Array1::zeros(steps),
        0.0001,
        100,
        n,
        e,
        tau,
    );
    calculate_f_from_eccentric_anomaly(eccentric_anomaly, e)
}

/// Calculates the true anomaly f for a 2 body problem from keplerian elements for a given total time split into a given number of steps.
///
/// One of the bodies lies at the origin and all the inputs are given with respect to this body, referred to as "the central body". The other body is referred to as "the rotating body"
///
/// **Inputs**:
///
/// elements: Keplerian elements of the 2-body system
///
/// mu: The gravitational parameter of the system. See [μ](`calculate_mu`).
///
/// total_time: The time to be simulated.
///
/// steps: The number of intervals the original interval [0, total_time] is split into.
///
/// **Output**: An array of true anomalies.
pub fn calculate_f_from_keplerian_elements(
    elements: &KeplerianElements,
    mu: f64,
    start_time: f64,
    end_time: f64,
    steps: usize,
) -> Array1<f64> {
    let n: f64 = calculate_n(mu, elements.a);
    let t: Array1<f64> = Array1::linspace(start_time, end_time, steps);
    let eccentric_anomaly: Array1<f64> = calculate_eccentric_anomaly_iteratively(
        t.clone(),
        Array1::zeros(steps),
        0.0001,
        100,
        n,
        elements.e,
        elements.tau,
    );
    calculate_f_from_eccentric_anomaly(eccentric_anomaly, elements.e)
}

/// Calculates the semi-major axis a for a 2 body problem from initial conditions.
///
/// **Inputs**:
///
/// rr: The initial position of the rotating body with respect to the central body.
///
/// vv: The initial velocity of the rotating body with respect to the central body.
///
/// mu: The gravitational parameter of the system. See [μ](`calculate_mu`).
///
/// **Output**: Semi-major axis a.
pub fn calculate_a_from_initial_rr_and_vv(rr: Array1<f64>, vv: Array1<f64>, mu: f64) -> f64 {
    let h: f64 = calculate_h(rr, vv, mu);
    calculate_a(mu, h)
}
