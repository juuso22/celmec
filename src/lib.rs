/// This module contains functions to calculate various properties of orbits.
///
/// It is currently limited to deal with the 2-body problem.
///
/// ## Naming
/// Naming tends to cause headaches at times and this module is no different. Here are some terms I use in a way that are to my best knowledge non-standard, but for which the standard terminology escapes me:
///
/// Eccentric anomaly: I call the equivalent of the eccentric anomaly for elliptic orbit eccentric anomaly also in the context of prabolic and hyperbolic orbits.
///
/// 'hyperbolic Kepler's equation': This is the equation M = E - e * sinh(E) which relates mean anomaly M to eccentric anomaly (see note above) for hyperbolic orbits (with help from eccentricity e). I donät know what the correct name for this equation is.
pub mod orbit;

/// This module contains structs for orbital elements.
///
/// It is currently limited to deal with the 2-body problem.
pub mod orbital_elements;

/// This module contains functions for transformations between coordinate systems.
///
/// Eg. switching from the orbital plane of a 2-body system to spherical coordinates with one of the bodies at origin.
pub mod transformations;

#[cfg(test)]
mod math_tests {
    use super::*;
    use ndarray::{array, Array1};
    use orbit::math::*;
    use std::collections::HashMap;

    #[test]
    fn non_converging_iterations() {
        fn square(init_v: Array1<f64>, t: Array1<f64>, _params: HashMap<&str, f64>) -> Array1<f64> {
            init_v.mapv_into(|v| v.powf(2.)) + t
        }
        assert_eq!(
            solve_equation_iteratively(
                &square,
                array![2., 3.],
                array![0., 0.],
                HashMap::new(),
                0.1,
                4
            ),
            array![65536., 43046721.]
        )
    }

    #[test]
    fn converging_iterations() {
        fn square(init_v: Array1<f64>, t: Array1<f64>, _params: HashMap<&str, f64>) -> Array1<f64> {
            -init_v.mapv_into(|v| v.powf(2.)) + t + 1.
        }
        let tolerance: f64 = 0.00000001;
        assert!(
            solve_equation_iteratively(
                &square,
                array![0.1],
                array![0.],
                HashMap::new(),
                tolerance,
                99
            )[0] - 1.
                < tolerance
        )
    }

    #[test]
    fn cross_and_dot_product_combine_to_almost_zero() {
        let a: Array1<f64> = array![2., -1.2, 1.];
        let b: Array1<f64> = array![-5., 0.00023, -3.];
        assert!((cross_product(a.clone(), b) * a).sum() < 1.0e-10)
    }

    #[test]
    fn eculidean_norm_of_three_and_four_is_five() {
        assert_eq!(euclidean_norm(array![3., 4., 0.]), 5.)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{array, Array1};
    use std::f64::consts::PI;

    #[test]
    fn k_dot_e_is_almost_zero() {
        let rr0: Array1<f64> = array![1., 3.2, 0.4];
        let v0: Array1<f64> = array![-4., 0.002, 2.2];
        let mu: f64 = 5.123;
        let ee: Array1<f64> = orbit::calculate_ee(rr0.clone(), v0.clone(), mu);
        let kk: Array1<f64> = orbit::math::cross_product(rr0, v0);
        assert!((ee * kk).sum() < 1.0e-10)
    }

    #[test]
    fn earth_e() {
        let solar_mass: f64 = 1.98847 * 10_f64.powf(30.);
        let earth_mass: f64 = 5.972 * 10_f64.powf(24.);
        let mu: f64 = orbit::calculate_mu(solar_mass, earth_mass);
        let earth_perihelion = array![147.10 * 10_f64.powf(9.), 0., 0.];
        let earth_orbital_vv_at_perihelion = array![0., 30290., 0.];
        assert!(
            orbit::calculate_e(
                earth_perihelion.clone(),
                earth_orbital_vv_at_perihelion.clone(),
                mu
            ) - 0.0167
                < 0.0005
        );
        assert!(
            orbit::calculate_e(earth_perihelion, earth_orbital_vv_at_perihelion, mu) - 0.0167
                > -0.0005
        );
    }

    #[test]
    fn simple_mean_anomaly() {
        assert_eq!(
            array![-4., -5.],
            orbit::calculate_mean_anomaly(array![1., 0.5], 2., 3.)
        );
    }

    #[test]
    fn simple_f_from_series() {
        assert_eq!(
            array![PI, PI / 2. + 1.75 - 13. / 12.],
            orbit::calculate_f_from_series(array![1., 0.5], 1., 2., 0.)
        );
    }

    #[test]
    fn radius_from_f_with_zero_e() {
        let a: f64 = 1.;
        assert_eq!(array![a], orbit::calculate_radius_from_f(array![3.], 0., a));
    }

    #[test]
    fn polar_angle_for_zero_iota() {
        let iota: f64 = 0.;
        let omega: f64 = 0.;
        assert_eq!(
            array![PI / 2., PI / 2., PI / 2.],
            transformations::polar_angle_from_keplerian_elements(array![1., 2., 3.], iota, omega)
        );
    }

    #[test]
    fn azimuthal_angle_for_zero_iota() {
        let iota: f64 = 0.;
        let omega: f64 = 0.;
        let longitude_of_the_ascending_node: f64 = 0.;
        let f: Array1<f64> = array![0., PI / 2., PI, 3. / 2. * PI];
        let polar_angle: Array1<f64> =
            transformations::polar_angle_from_keplerian_elements(f.clone(), iota, omega);
        assert_eq!(
            f.clone(),
            transformations::azimuthal_angle_from_keplerian_elements_and_polar_angle(
                f,
                polar_angle,
                omega,
                longitude_of_the_ascending_node,
            )
        );
    }

    #[test]
    fn polar_angle_for_pi_per_two_iota() {
        let iota: f64 = PI / 2.;
        let omega: f64 = 0.;
        let f = orbit::calculate_f_from_series(Array1::linspace(0., 1., 10), 0., 1., 0.);
        let f_sin = f.clone().mapv_into(|v| v.sin());
        let polar_angle_opposite =
            PI / 2. - transformations::polar_angle_from_keplerian_elements(f.clone(), iota, omega);
        let polar_angle_opposite_sin = polar_angle_opposite.clone().mapv_into(|v| v.sin());
        assert!((f_sin - polar_angle_opposite_sin)
            .mapv_into_any(|v| v < 0.00000001 && v > -0.00000001)
            .fold(true, |a, b| a && *b));
        let f_cos = f.clone().mapv_into(|v| v.cos());
        let polar_angle_opposite_cos = polar_angle_opposite.clone().mapv_into(|v| v.cos());
        assert!((f_cos.mapv_into(|v| v.abs()) - polar_angle_opposite_cos)
            .mapv_into_any(|v| v < 0.00000001 && v > -0.00000001)
            .fold(true, |a, b| a && *b));
    }

    #[test]
    fn azimuthal_angle_for_pi_longitude_of_ascending_node() {
        let omega: f64 = 0.;
        let iota: f64 = PI / 4.;
        let f: Array1<f64> =
            orbit::calculate_f_from_series(Array1::linspace(0., 1., 10), 0., 1., 0.);
        let polar_angle: Array1<f64> =
            transformations::polar_angle_from_keplerian_elements(f.clone(), iota, omega);
        let azimuthal_angle_zero_loan: Array1<f64> =
            transformations::azimuthal_angle_from_keplerian_elements_and_polar_angle(
                f.clone(),
                polar_angle.clone(),
                omega,
                0.,
            );
        let azimuthal_angle_pi_loan: Array1<f64> =
            transformations::azimuthal_angle_from_keplerian_elements_and_polar_angle(
                f.clone(),
                polar_angle.clone(),
                omega,
                PI,
            );
        assert_eq!(azimuthal_angle_zero_loan + PI, azimuthal_angle_pi_loan);
    }

    #[test]
    fn azimuthal_angle_for_pi_omega() {
        let longitude_of_ascending_node: f64 = 0.;
        let iota: f64 = PI / 4.;
        let omega: f64 = 0.;
        let half_ticks: usize = 5;
        let ticks: usize = 2 * half_ticks;
        let f: Array1<f64> =
            orbit::calculate_f_from_series(Array1::linspace(0., 1., ticks), 0., 1., 0.);
        let polar_angle: Array1<f64> =
            transformations::polar_angle_from_keplerian_elements(f.clone(), iota, omega);
        let azimuthal_angle_zero_aop: Array1<f64> =
            transformations::azimuthal_angle_from_keplerian_elements_and_polar_angle(
                f.clone(),
                polar_angle.clone(),
                0.,
                longitude_of_ascending_node,
            );
        let azimuthal_angle_pi_aop: Array1<f64> =
            transformations::azimuthal_angle_from_keplerian_elements_and_polar_angle(
                f.clone(),
                polar_angle.clone(),
                PI,
                longitude_of_ascending_node,
            );
        let iter1 = azimuthal_angle_pi_aop.iter();
        for (i, v) in iter1.enumerate() {
            println!("{}", (i + half_ticks) % ticks);
            assert_eq!(*v, azimuthal_angle_zero_aop[(i + half_ticks) % ticks]);
        }
        assert!(azimuthal_angle_pi_aop
            .iter()
            .enumerate()
            .map(|(idx, v)| *v == azimuthal_angle_zero_aop[(idx + half_ticks) % ticks])
            .fold(true, |a, b| a && b));
    }
}
