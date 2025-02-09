use crate::math::{atan, cross_product, euclidean_norm};
use crate::mechanics::calculate_kk_from_initial_rr_and_vv;
use crate::two_body::{
    calculate_eccentric_anomaly_from_f, calculate_h, calculate_initial_f_from_initial_conditions,
    calculate_mean_anomaly_from_eccentric_anomaly, calculate_n,
};
use ndarray::{array, Array1};

/// Keplerian elements of a 2-body system
///
/// **Fields**:
///
/// e = eccentricity
///
/// longitude_of_the_ascending_node = longitude of the ascending node
///
/// tau = perihelion time
///
/// a = semi-major axis
///
/// i = inclination
///
/// omega = argument of perihelion
#[derive(Debug, Clone)]
pub struct KeplerianElements {
    pub e: f64,
    pub a: f64,
    pub tau: f64,
    pub iota: f64,
    pub omega: f64,
    pub longitude_of_the_ascending_node: f64,
}

/// Calculates inclination from initial position and velocity
///
/// **Inputs**
///
/// rr: The initial position of the rotating body with respect to the central body.
///
/// vv: The initial velocity of the rotating body with respect to the central body.
///
/// **Output**: inclination of an orbit
pub fn calculate_iota_from_initial_rr_and_vv(rr: Array1<f64>, vv: Array1<f64>) -> f64 {
    let kk: Array1<f64> = calculate_kk_from_initial_rr_and_vv(rr.clone(), vv.clone());
    calculate_iota_from_kk(kk)
}

/// Calculates inclination from angular momentum per unit mass
///
/// **Inputs**
///
/// kk: angular momentum per unit mass
///
/// **Output**: inclination of an orbit
pub fn calculate_iota_from_kk(kk: Array1<f64>) -> f64 {
    let k: f64 = euclidean_norm(kk.clone());
    (kk[2] / k).acos()
}

/// Calculates longitude of the ascending node from initial position and velocity
///
/// **Inputs**
///
/// rr: The initial position of the rotating body with respect to the central body.
///
/// vv: The initial velocity of the rotating body with respect to the central body.
///
/// **Output**: longitude of the ascending node of an orbit
pub fn calculate_longitude_of_the_ascending_node_from_initial_rr_and_vv(
    rr: Array1<f64>,
    vv: Array1<f64>,
) -> f64 {
    let kk: Array1<f64> = calculate_kk_from_initial_rr_and_vv(rr, vv);
    calculate_longitude_of_the_ascending_node_from_kk(kk)
}

/// Calculates longitude of the ascending node from angular momentum per unit mass vector and inclination
///
/// **Inputs**
///
/// kk: anuglar momentum per unit mass of the system.
///
/// **Output**: longitude of the ascending node of an orbit
pub fn calculate_longitude_of_the_ascending_node_from_kk(kk: Array1<f64>) -> f64 {
    let ascending_node_vector: Array1<f64> =
        kk[2].signum() * cross_product(array![0., 0., kk[2]], kk);
    atan(ascending_node_vector[1], ascending_node_vector[0])
}

/// Calculates argument of perihelion (omega) from the gravitational paramerter and some initial position and velocity
///
/// **Inputs**
///
/// mu: gravitational parameter. See [calculate_mu](`calculate_mu`).
///
/// rr: The initial position of the rotating body with respect to the central body.
///
/// vv: The initial velocity of the rotating body with respect to the central body.
///
/// **Output**: argument of perihelion (omega) of an orbit
pub fn calculate_omega_from_mu_and_initial_rr_and_vv(
    mu: f64,
    rr: Array1<f64>,
    vv: Array1<f64>,
) -> f64 {
    let kk: Array1<f64> = calculate_kk_from_initial_rr_and_vv(rr.clone(), vv.clone());
    let ee: Array1<f64> = calculate_ee(rr, vv, mu);
    calculate_omega_from_kk_and_ee(kk, ee)
}

/// Calculates argument of perihelion (omega) from angular momentum per unit mass vector, the "eccentricity vector", eccentrcity and inclination
///
/// **Inputs**
///
/// kk: anuglar momentum per unit mass of the system.
///
/// ee: eccentricity_vector
///
/// **Output**: argument of perihelion (omega) of an orbit
pub fn calculate_omega_from_kk_and_ee(kk: Array1<f64>, ee: Array1<f64>) -> f64 {
    if (kk[0] != 0.) || (kk[1] != 0.) {
        let ascending_node_vector: Array1<f64> =
            kk[2].signum() * cross_product(array![0., 0., kk[2]], kk);
        let cos_omega: f64 = (ascending_node_vector.clone() * ee.clone())
            .iter()
            .sum::<f64>();
        let sin_omega: f64 = euclidean_norm(cross_product(ascending_node_vector, ee));
        atan(sin_omega, cos_omega)
    } else {
        atan(ee[1], ee[0])
    }
}

/// Calculates vector **e** between 2 bodies.
///
/// This is the constant vector whose length gives the [eccentricity e](`calculate_e`) of an orbit for the 2-body problem.
///
/// **e** = - (**r** x **v**) / μ - **r** / r,
///
/// where
///
/// **r** = position of the 2 bodies relative to each other at some point in time,
///
/// **v** = velocity of the 2 bodies relative to each other at the same point of time as **r**
///
/// μ = see [μ](`calculate_mu`)
///
/// r = |**r**|
///
/// Inputs are a position **r** (`rr`) and the velocity **v** (`vv`) of the 2 bodies with respect to each other at any one point of time as well as their [μ](`calculate_mu`).
pub fn calculate_ee(rr: Array1<f64>, vv: Array1<f64>, mu: f64) -> Array1<f64> {
    -cross_product(
        calculate_kk_from_initial_rr_and_vv(rr.clone(), vv.clone()),
        vv,
    ) / mu
        - rr.clone() / euclidean_norm(rr)
}

/// Calculates eccentricity of an orbit for 2 bodies.
///
/// **Inputs**:
///
/// rr = position of the 2 bodies with respect to each other
///
/// vv = velocity of the 2 bodies with respect to each other
///
/// μ = see [μ](`calculate_mu`).
///
/// **Output**: eccentricity e
///
/// Calculation is done from e's definition as the length of the vector **e** ([`ee`](`calculate_ee`)).
pub fn calculate_e(rr: Array1<f64>, vv: Array1<f64>, mu: f64) -> f64 {
    euclidean_norm(calculate_ee(rr, vv, mu))
}

/// Calculates the semi-major axis a of a 2-body system
///
/// a = μ / (2 * h), if h >= 0 and -μ / (2 * h) otherwise,
///
/// where
///
/// μ = see [μ](`calculate_mu`)
///
/// h = Lagrangian of the system ie. the difference between its kinetic and potential energy.
///
/// Inputs are [μ](`calculate_mu`) and Lagrangian [h](`calculate_h`) of the system.
pub fn calculate_a(mu: f64, h: f64) -> f64 {
    let mut a: f64 = mu / 2. / h;
    if h < 0. {
        a = -1. * a;
    }
    a
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
pub fn calculate_tau_from_mean_anomaly(t: f64, mean_anomaly: f64, n: f64) -> f64 {
    t - mean_anomaly / n
}

/// Calculates perihelion time of a 2-body system from initial conditions
///
/// The calculation is done so that first initial true anomaly f is calculated from the initial position and velocity. F is then used to calculate initial eccentric anomaly, which in turn is used to calculate initial mean anomaly. From the mean anomaly, tau is determined using the time when the initial position and velocity occur.
///
///
/// **Inputs**:
///
/// rr: position of the 2 bodies with respect to each other
///
/// vv: velocity of the 2 bodies with respect to each other
///
/// mu: see [μ](`calculate_mu`).
///
/// start_time: the time at which rr and vv occur
///
/// **Output**: Perihelion time of a 2-body system
pub fn calculate_tau_from_initial_rr_and_vv_mu_and_start_time(
    rr: Array1<f64>,
    vv: Array1<f64>,
    mu: f64,
    start_time: f64,
) -> f64 {
    let ee: Array1<f64> = calculate_ee(rr.clone(), vv.clone(), mu);
    let e: f64 = calculate_e(rr.clone(), vv.clone(), mu);
    let a: f64 = calculate_a_from_initial_rr_and_vv(rr.clone(), vv, mu);
    let initial_f: f64 = calculate_initial_f_from_initial_conditions(rr, ee, e);
    let initial_eccentric_anomaly: Array1<f64> =
        calculate_eccentric_anomaly_from_f(array![initial_f], e);
    let n: f64 = calculate_n(mu, a);
    calculate_tau_from_mean_anomaly(
        start_time,
        calculate_mean_anomaly_from_eccentric_anomaly(initial_eccentric_anomaly, e)[0],
        n,
    )
}

/// Calculates perihelion time of a 2-body system from initial conditions
///
/// The calculation is done so that first initial true anomaly f is calculated from the initial position and velocity. F is then used to calculate initial eccentric anomaly, which in turn is used to calculate initial mean anomaly. From the mean anomaly, tau is determined using the time when the initial position and velocity occur.
///
///
/// **Inputs**:
///
/// rr: position of the 2 bodies with respect to each other
///
/// mu: see [μ](`calculate_mu`).
///
/// start_time: the time at which rr and vv occur
///
/// **Output**: Perihelion time of a 2-body system
pub fn calculate_tau_from_initial_rr_mu_start_time_a_e_and_ee(
    rr: Array1<f64>,
    mu: f64,
    start_time: f64,
    a: f64,
    e: f64,
    ee: Array1<f64>,
) -> f64 {
    let initial_f: f64 = calculate_initial_f_from_initial_conditions(rr, ee.clone(), e);
    let initial_eccentric_anomaly: Array1<f64> =
        calculate_eccentric_anomaly_from_f(array![initial_f], e);
    let n: f64 = calculate_n(mu, a);
    calculate_tau_from_mean_anomaly(
        start_time,
        calculate_mean_anomaly_from_eccentric_anomaly(initial_eccentric_anomaly, e)[0],
        n,
    )
}

/// Calculates Keplerien elements other than perihelion time from an initial position and velocity and μ of a 2-body system and the time at which the initial conditions occur
///
/// **Inputs**:
///
/// rr: position of the 2 bodies with respect to each other
///
/// vv: velocity of the 2 bodies with respect to each other
///
/// mu: see [μ](`calculate_mu`).
///
/// start_time: the time at which rr and vv occur
///
/// **Output**: Keplerian elements of a 2-body system
pub fn calculate_keplerian_elements_from_initial_rr_and_vv_and_mu(
    rr: Array1<f64>,
    vv: Array1<f64>,
    mu: f64,
    start_time: f64,
) -> KeplerianElements {
    let kk: Array1<f64> = calculate_kk_from_initial_rr_and_vv(rr.clone(), vv.clone());
    let ee: Array1<f64> = calculate_ee(rr.clone(), vv.clone(), mu);
    let e: f64 = calculate_e(rr.clone(), vv.clone(), mu);
    let h: f64 = calculate_h(rr.clone(), vv.clone(), mu);
    let a: f64 = calculate_a(mu, h);
    KeplerianElements {
        e,
        longitude_of_the_ascending_node: calculate_longitude_of_the_ascending_node_from_kk(
            kk.clone(),
        ),
        tau: calculate_tau_from_initial_rr_mu_start_time_a_e_and_ee(
            rr,
            mu,
            start_time,
            a,
            e,
            ee.clone(),
        ),
        a,
        iota: calculate_iota_from_kk(kk.clone()),
        omega: calculate_omega_from_kk_and_ee(kk, ee),
    }
}
