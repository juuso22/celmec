use ndarray::Array1;
use std::f64::consts::PI;

pub fn polar_angle_from_keplerian_elements(
    f: Array1<f64>,
    iota: f64,
    omega: f64,
) -> Array1<f64> {
    PI / 2.
        - ((f + omega).mapv_into(|f| f.sin()) * iota.sin())
            .mapv_into(|f| f.asin())
}

fn azimuthal_angle_refinment(
    raw_azimuthal_angle: Array1<f64>,
    polar_angle: Array1<f64>,
) -> Array1<f64> {
    let mut azimuthal_angle: Array1<f64> = raw_azimuthal_angle.clone();
    for (i, v) in raw_azimuthal_angle.iter().enumerate() {
        if polar_angle[i] >= PI / 2. {
            azimuthal_angle[i] = -*v;
        }
    }
    azimuthal_angle
}

pub fn azimuthal_angle_from_keplerian_elements_and_polar_angle(
    f: Array1<f64>,
    polar_angle: Array1<f64>,
    omega: f64,
    longitude_of_the_ascending_node: f64,
) -> Array1<f64> {
    let angle_from_ascending_node = f.clone() + omega;
    let mut polar_angle_iterator = polar_angle.iter();
    let raw_azimuthal_angle: Array1<f64> = angle_from_ascending_node
        .mapv_into(|f| f.cos() / ((PI / 2.) - polar_angle_iterator.next().unwrap()).cos())
        .mapv_into(|f| f.acos());
    azimuthal_angle_refinment(raw_azimuthal_angle, polar_angle) + longitude_of_the_ascending_node
}
