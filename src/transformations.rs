use ndarray::Array1;
use std::f64::consts::PI;

/// Calculates polar angle from true anomaly f and relevant Keplerian elements.
///
/// **Inputs**:
///
/// f: array of true anomalies
///
/// iota: inclination
///
/// omega: argument of perihelion
///
/// **Output**: array of polar angles
///
/// Polar angle is calculated in a reference frame where one of the two bodies of a 2-body system is at orign from the formula:
///
/// θ = π / 2 - asin(sin(f + ω) * sin(ι)),
///
/// where
///
/// θ = polar angle
///
/// f = true anomaly
///
/// ω = argument of perihelion
///
/// ι = inclination
pub fn theta_from_keplerian_elements(f: Array1<f64>, iota: f64, omega: f64) -> Array1<f64> {
    PI / 2. - ((f + omega).mapv_into(|f| f.sin()) * iota.sin()).mapv_into(|f| f.asin())
}

/// Determine whether the azimuthal angle is in range [-π, 0] or (0, π] based on polar angle.
///
/// **Inputs**
///
/// raw_phi: array azimuthal angles in range [0, π]
///
/// polar_angel: array of polar angle corrensponding to the "raw" aimuthal angles
///
/// **Output**: an array of azimuthal angles in range [-π, π]
///
/// An azimuthal angle is negated if the corresponding polar angle is greater or equal than π/2.
fn phi_refinment(raw_phi: Array1<f64>, theta: Array1<f64>) -> Array1<f64> {
    let mut phi: Array1<f64> = raw_phi.clone();
    for (i, v) in raw_phi.iter().enumerate() {
        if theta[i] >= PI / 2. {
            phi[i] = -*v;
        }
    }
    phi
}

/// Calculates azimuthal angle from known polar angles, true anomalies, and relevant Keplerian elements (argument of perihelion and longitude of the ascending node).
///
/// **Inputs**:
///
/// f: an array of true anomalies
///
/// theta: an array of polar angles corresponfing to the true anomalies f
///
/// omega: argument of perihelion
///
/// longitude_of_the_ascending_node: longitude of the ascending node
///
/// **Output**: an array of azimuthal angles corresponding to the true anomalies f (and polar angles theta)
///
/// Following formula is used for the calculation:
///
/// φ = acos(cos(f + ω) / (π / 2) - cos(θ)) + Ω,
///
/// where
///
/// φ = azimutahl angle
///
/// θ = polar angle
///
/// f = true anomaly
///
/// ω = argument of perihelion
///
/// Ω = longitude of the ascending node
///
/// As acos returns angles in range [0, π], polar angle is used to flip the azimuthal angle's sign if needed (ie. if θ >= π/2)
pub fn phi_from_keplerian_elements_and_theta(
    f: Array1<f64>,
    theta: Array1<f64>,
    omega: f64,
    longitude_of_the_ascending_node: f64,
) -> Array1<f64> {
    let angle_from_ascending_node = f.clone() + omega;
    let mut theta_iterator = theta.iter();
    let raw_phi: Array1<f64> = angle_from_ascending_node
        .mapv_into(|f| f.cos() / ((PI / 2.) - theta_iterator.next().unwrap()).cos())
        .mapv_into(|f| f.acos());
    phi_refinment(raw_phi, theta) + longitude_of_the_ascending_node
}
