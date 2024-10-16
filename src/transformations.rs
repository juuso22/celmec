use ndarray::Array1;
use std::f64::consts::PI;

use crate::orbital_elements;

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
/// Polar angle is calculated in a reference frame where one of the two bodies of a 2-body system is at origin from the formula:
///
/// θ = atan(sin(f + ω) * cos(ι) / cos(f + ω)) + Ω,
///
/// where
///
/// θ = polar angle
///
/// f = true anomaly
///
/// ι = inclination
///
/// ω = argument of perihelion
///
/// Ω = longitude of the ascending node
///
/// Because asin values are bound between -π/2 and π/2, the resulting polar angles need to be shifted to the second and third qudrant of the unit circle if needed. The value of of f+ω is always is the same quadrant of the unit circle as value projected to the epliptical plane (ie. the value that comes out of the asin, before adding the longitude of the ascending node in the formula above). Therefore, if that value is in the second or third quadrant, the value coming out of the asin is shifted from first quadrant to second or fourth quadrant to the third:
///
/// θ<sub>projected</sub> = atan(sin(f + ω) * cos(ι) / cos(f + ω)), if f+ω ∈ [-π/2, π/2]
///
/// θ<sub>projected</sub> = π + atan(sin(f + ω) * cos(ι) / cos(f + ω)), otherwise
///
/// where
///
/// θ<sub>projected</sub> = the polar angle projected on the ecliptica before adding longitude of the ascending node
///
/// In the end, the polar angle is bound to the interval [0, 2π).
pub fn theta_from_keplerian_elements(
    f: Array1<f64>,
    iota: f64,
    omega: f64,
    longitude_of_the_ascending_node: f64,
) -> Array1<f64> {
    let f_plus_omega_bound: Array1<f64> =
        ((f + omega) % (2. * PI)).mapv_into(|v| if v < 0. { 2. * PI + v } else { v });
    let over_pi_per_two_mask: Array1<f64> = (f_plus_omega_bound.clone())
        .mapv_into(|f| ((f > PI / 2.) & (f < PI * 3. / 2.)) as i32 as f64);
    ((PI * over_pi_per_two_mask
        + (f_plus_omega_bound).mapv_into_any(|f| {
            if f.sin() == -1. {
                let iota_norm = iota % (2. * PI);
                if (iota_norm == (PI / 2.)) || (iota_norm == (3. * PI / 2.)) {
                    f64::NAN
                } else {
                    3. * PI / 2.
                }
            } else if f.sin() == 1. {
                let iota_norm = iota % (2. * PI);
                if (iota_norm == (PI / 2.)) || (iota_norm == (3. * PI / 2.)) {
                    f64::NAN
                } else {
                    PI / 2.
                }
            } else {
                ((f.sin() * iota.cos().abs()) / f.cos()).atan()
            }
        })
        + longitude_of_the_ascending_node)
        % (2. * PI))
        .mapv_into(|v| if v < 0. { 2. * PI + v } else { v })
}

/// Calculates azimuthal angle from known true anomalies and relevant Keplerian elements (inclination and argument of perihelion).
///
/// **Inputs**:
///
/// f: an array of true anomalies
///
/// iota: inclination
///
/// omega: argument of perihelion
///
/// **Output**: an array of azimuthal angles corresponding to the true anomalies f
///
/// Following formula is used for the calculation:
///
/// φ = asin(sin(f + ω) * sin(ι)),
///
/// where
///
/// φ = azimutahl angle
///
/// f = true anomaly
///
/// ι = inclination
///
/// ω = argument of perihelion
pub fn phi_from_keplerian_elements(f: Array1<f64>, iota: f64, omega: f64) -> Array1<f64> {
    PI / 2. - ((f + omega).mapv_into(|v| v.sin()) * iota.sin()).mapv_into(|v| v.asin())
}

/// Calculates spherical angle coordinates theta and phi corresponding to ecliptical coordinates from an array of true anomalies and orbital elements
///
/// **Inputs:***
///
/// f: an array of true anomalies
///
/// keplerian_elements: a struct containing Keplerian orbital elements
///
/// **Output:** A tuple where the first array is an array of sphrical angles theta and the the second array is an array of spherical angles phi
///
/// The math behind this method can be found from the functions it uses: [theta_from_keplerian_elements](`theta_from_keplerian_elements`) and [phi_from_keplerian_elements](`phi_from_keplerian_elements`)
pub fn spherical_coordinates_from_f_and_keplerian_elements(
    f: Array1<f64>,
    keplerian_elements: orbital_elements::KeplerianElements,
) -> (Array1<f64>, Array1<f64>) {
    let theta: Array1<f64> = theta_from_keplerian_elements(
        f.clone(),
        keplerian_elements.iota,
        keplerian_elements.omega,
        keplerian_elements.longitude_of_the_ascending_node,
    );
    let phi: Array1<f64> =
        phi_from_keplerian_elements(f, keplerian_elements.iota, keplerian_elements.omega);
    (theta, phi)
}

/// Calculates cartesian coordinates from spherical coordinates
///
/// **Inputs:***
///
/// theta: an array of spherical angles theta
///
/// phi: an array of spherical angles phi
///
/// r: an array of radii
///
/// **Output:** A triple where the elements ar arrays of x, y and z coordinates, respectively
pub fn cartesian_coordinates_from_spherical_coordinates(
    theta: Array1<f64>,
    phi: Array1<f64>,
    r: Array1<f64>,
) -> (Array1<f64>, Array1<f64>, Array1<f64>) {
    let x: Array1<f64> =
        r.clone() * phi.clone().mapv_into(|v| v.sin()) * theta.clone().mapv_into(|v| v.cos());
    let y: Array1<f64> =
        r.clone() * phi.clone().mapv_into(|v| v.sin()) * theta.mapv_into(|v| v.sin());
    let z: Array1<f64> = r * phi.mapv_into(|v| v.cos());
    (x, y, z)
}

/// Calculates cartesian coordinates corresponding to the ecliptical coordinates from an array of true anomalies and orbital elements
///
/// **Inputs:***
///
/// f: an array of true anomalies
///
/// keplerian_elements: a struct containing Keplerian orbital elements
///
/// **Output:** A triple where the elements ar arrays of x, y and z coordinates, respectively
pub fn cartesian_coordinates_from_f_r_and_keplerian_elements(
    f: Array1<f64>,
    r: Array1<f64>,
    keplerian_elements: orbital_elements::KeplerianElements,
) -> (Array1<f64>, Array1<f64>, Array1<f64>) {
    let (theta, phi) = spherical_coordinates_from_f_and_keplerian_elements(f, keplerian_elements);
    cartesian_coordinates_from_spherical_coordinates(theta, phi, r)
}
