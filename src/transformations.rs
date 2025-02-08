use ndarray::{Array1, Array2};
use std::f64::consts::PI;

use crate::orbital_elements;

/// Calculates azimuthal angle from true anomaly f and relevant Keplerian elements.
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
/// φ = atan(sin(f + ω) * cos(ι) / cos(f + ω)) + Ω,
///
/// where
///
/// φ = azimuthal angle
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
/// φ<sub>projected</sub> = atan(sin(f + ω) * cos(ι) / cos(f + ω)), if f+ω ∈ [-π/2, π/2]
///
/// φ<sub>projected</sub> = π + atan(sin(f + ω) * cos(ι) / cos(f + ω)), otherwise
///
/// where
///
/// φ<sub>projected</sub> = the polar angle projected on the ecliptica before adding longitude of the ascending node
///
/// In the end, the polar angle is bound to the interval [0, 2π).
pub fn phi_from_keplerian_elements(
    f: Array1<f64>,
    iota: f64,
    omega: f64,
    longitude_of_the_ascending_node: f64,
) -> Array1<f64> {
    let mut non_nan_longitude_of_the_ascending_node: f64 = longitude_of_the_ascending_node;
    if non_nan_longitude_of_the_ascending_node.is_nan() {
        non_nan_longitude_of_the_ascending_node = 0.;
    }
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
        + non_nan_longitude_of_the_ascending_node)
        % (2. * PI))
        .mapv_into(|v| if v < 0. { 2. * PI + v } else { v })
}

/// Calculates polar angle from known true anomalies and relevant Keplerian elements (inclination and argument of perihelion).
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
/// θ = asin(sin(f + ω) * sin(ι)),
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
pub fn theta_from_keplerian_elements(f: Array1<f64>, iota: f64, omega: f64) -> Array1<f64> {
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
/// **Output:** An array of arrays where the first array is an array of sphrical angles theta and the the second array is an array of spherical angles phi
///
/// The math behind this method can be found from the functions it uses: [theta_from_keplerian_elements](`theta_from_keplerian_elements`) and [phi_from_keplerian_elements](`phi_from_keplerian_elements`)
pub fn spherical_coordinates_from_f_and_keplerian_elements(
    f: Array1<f64>,
    keplerian_elements: orbital_elements::KeplerianElements,
) -> Array2<f64> {
    let phi: Array1<f64> = phi_from_keplerian_elements(
        f.clone(),
        keplerian_elements.iota,
        keplerian_elements.omega,
        keplerian_elements.longitude_of_the_ascending_node,
    );
    let theta: Array1<f64> =
        theta_from_keplerian_elements(f, keplerian_elements.iota, keplerian_elements.omega);
    let mut spherical_coordinates: Array2<f64> = Array2::zeros((2, theta.len()));
    spherical_coordinates.row_mut(0).assign(&phi);
    spherical_coordinates.row_mut(1).assign(&theta);
    spherical_coordinates
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
/// **Output:** An array of arrays where the inner arrays are x, y and z coordinates, respectively
pub fn cartesian_coordinates_from_spherical_coordinates(
    phi: Array1<f64>,
    theta: Array1<f64>,
    r: Array1<f64>,
) -> Array2<f64> {
    let mut xyz: Array2<f64> = Array2::zeros((3, r.len()));
    xyz.row_mut(0).assign(
        &(r.clone() * theta.clone().mapv_into(|v| v.sin()) * phi.clone().mapv_into(|v| v.cos())),
    );
    xyz.row_mut(1)
        .assign(&(r.clone() * theta.clone().mapv_into(|v| v.sin()) * phi.mapv_into(|v| v.sin())));
    xyz.row_mut(2).assign(&(r * theta.mapv_into(|v| v.cos())));
    xyz
}

/// Calculates cartesian coordinates from polar coordinate.
///
/// **Inputs:***
///
/// theta: an array of polar angles theta
///
/// r: an array of radii
///
/// **Output:** An array of arrays where the inner arrays are x, y and z coordinates, respectively
///
/// Z is assumed to be 0.
pub fn cartesian_coordinates_from_polar_coordinates(
    theta: Array1<f64>,
    r: Array1<f64>,
) -> Array2<f64> {
    let mut xyz: Array2<f64> = Array2::zeros((3, r.len()));
    xyz.row_mut(0)
        .assign(&(r.clone() * theta.clone().mapv_into(|v| v.cos())));
    xyz.row_mut(1).assign(&(r * theta.mapv_into(|v| v.sin())));
    xyz
}

/// Calculates cartesian coordinates corresponding to the ecliptical coordinates from an array of true anomalies and orbital elements
///
/// **Inputs:***
///
/// f: an array of true anomalies
///
/// keplerian_elements: a struct containing Keplerian orbital elements
///
/// **Output:** An array of arrays where the inner arrays are x, y and z coordinates, respectively
pub fn cartesian_coordinates_from_f_r_and_keplerian_elements(
    f: Array1<f64>,
    r: Array1<f64>,
    keplerian_elements: orbital_elements::KeplerianElements,
) -> Array2<f64> {
    if keplerian_elements.longitude_of_the_ascending_node.is_nan() {
        let polar_angle: Array1<f64> = f - keplerian_elements.omega;
        cartesian_coordinates_from_polar_coordinates(polar_angle, r)
    } else {
        let spherical_coordinates: Array2<f64> =
            spherical_coordinates_from_f_and_keplerian_elements(f, keplerian_elements);
        cartesian_coordinates_from_spherical_coordinates(
            spherical_coordinates.row(0).to_owned(),
            spherical_coordinates.row(1).to_owned(),
            r,
        )
    }
}
