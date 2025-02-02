use crate::math::{atan, cross_product, euclidean_norm};
use crate::mechanics::calculate_kk_from_initial_rr_and_vv;
use crate::orbit::calculate_h;
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
#[derive(Clone)]
pub struct KeplerianElements {
    pub e: f64,
    pub longitude_of_the_ascending_node: f64,
    pub tau: f64,
    pub a: f64,
    pub iota: f64,
    pub omega: f64,
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
    let ascending_node_vector: Array1<f64> =
        kk[2].signum() * cross_product(array![0., 0., kk[2]], kk);
    let cos_omega: f64 = (ascending_node_vector.clone() * ee.clone())
        .iter()
        .sum::<f64>();
    let sin_omega: f64 = euclidean_norm(cross_product(ascending_node_vector, ee));
    atan(sin_omega, cos_omega)
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

/// Calculates Keplerien elements other than perihelion time from an initial position and velocity and μ of a 2-body system
///
/// **Inputs**:
///
/// rr = position of the 2 bodies with respect to each other
///
/// vv = velocity of the 2 bodies with respect to each other
///
/// μ = see [μ](`calculate_mu`).
///
/// **Output**: Keplerian elements of a 2-body system
///
/// Perihelion time is set to 0. This function will eventually be modified so that perihelion time is calculated by setting the time of the initial **r** and **v** to zero.
pub fn calculate_keplerian_elements_from_initial_rr_and_vv_and_mu(
    rr: Array1<f64>,
    vv: Array1<f64>,
    mu: f64,
) -> KeplerianElements {
    let kk: Array1<f64> = calculate_kk_from_initial_rr_and_vv(rr.clone(), vv.clone());
    let ee: Array1<f64> = calculate_ee(rr.clone(), vv.clone(), mu);
    let e: f64 = calculate_e(rr.clone(), vv.clone(), mu);
    let h: f64 = calculate_h(rr.clone(), vv.clone(), mu);
    KeplerianElements {
        e,
        longitude_of_the_ascending_node: calculate_longitude_of_the_ascending_node_from_kk(
            kk.clone(),
        ),
        tau: 0.,
        a: calculate_a(mu, h),
        iota: calculate_iota_from_kk(kk.clone()),
        omega: calculate_omega_from_kk_and_ee(kk, ee),
    }
}
