pub mod orbit;
pub mod orbital_elements;
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
        let r0: Array1<f64> = array![1., 3.2, 0.4];
        let v0: Array1<f64> = array![-4., 0.002, 2.2];
        let mu: f64 = 5.123;
        let e: Array1<f64> = orbit::calculate_e(r0.clone(), v0.clone(), mu);
        let k: Array1<f64> = orbit::math::cross_product(r0, v0);
        assert!((e * k).sum() < 1.0e-10)
    }

    #[test]
    fn earth_eccentricity() {
        let solar_mass: f64 = 1.98847 * 10_f64.powf(30.);
        let earth_mass: f64 = 5.972 * 10_f64.powf(24.);
        let mu: f64 = orbit::calculate_mu(solar_mass, earth_mass);
        let earth_periapsis = array![147.10 * 10_f64.powf(9.), 0., 0.];
        let earth_orbital_velocity_at_periapsis = array![0., 30290., 0.];
        assert!(
            orbit::calculate_eccentricity(
                earth_periapsis.clone(),
                earth_orbital_velocity_at_periapsis.clone(),
                mu
            ) - 0.0167
                < 0.0005
        );
        assert!(
            orbit::calculate_eccentricity(earth_periapsis, earth_orbital_velocity_at_periapsis, mu)
                - 0.0167
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
    fn simple_true_anomaly_from_series() {
        assert_eq!(
            array![PI, PI / 2. + 1.75 - 13. / 12.],
            orbit::calculate_true_anomaly_from_series(array![1., 0.5], 1., 2., 0.)
        );
    }

    #[test]
    fn radius_from_true_anomaly_with_zero_eccentricity() {
        let a: f64 = 1.;
        assert_eq!(
            array![a],
            orbit::calculate_radius_from_true_anomaly(array![3.], 0., a)
        );
    }

    #[test]
    fn polar_angle_for_zero_inclination() {
        let inclination: f64 = 0.;
        let argument_of_periapsis: f64 = 0.;
        assert_eq!(
            array![PI / 2., PI / 2., PI / 2.],
            transformations::polar_angle_from_keplerian_elements(
                array![1., 2., 3.],
                inclination,
                argument_of_periapsis
            )
        );
    }

    #[test]
    fn azimuthal_angle_for_zero_inclination() {
        let inclination: f64 = 0.;
        let argument_of_periapsis: f64 = 0.;
        let longitude_of_the_ascending_node: f64 = 0.;
        let true_anomaly: Array1<f64> = array![0., PI / 2., PI, 3. / 2. * PI];
        let polar_angle: Array1<f64> = transformations::polar_angle_from_keplerian_elements(
            true_anomaly.clone(),
            inclination,
            argument_of_periapsis,
        );
        assert_eq!(
            true_anomaly.clone(),
            transformations::azimuthal_angle_from_keplerian_elements_and_polar_angle(
                true_anomaly,
                polar_angle,
                argument_of_periapsis,
                longitude_of_the_ascending_node,
            )
        );
    }

    #[test]
    fn polar_angle_for_pi_per_two_inclination() {
        let inclination: f64 = PI / 2.;
        let argument_of_periapsis: f64 = 0.;
        let true_anomaly =
            orbit::calculate_true_anomaly_from_series(Array1::linspace(0., 1., 10), 0., 1., 0.);
        let true_anomaly_sin = true_anomaly.clone().mapv_into(|v| v.sin());
        let polar_angle_opposite = PI / 2.
            - transformations::polar_angle_from_keplerian_elements(
                true_anomaly.clone(),
                inclination,
                argument_of_periapsis,
            );
        let polar_angle_opposite_sin = polar_angle_opposite.clone().mapv_into(|v| v.sin());
        assert!((true_anomaly_sin - polar_angle_opposite_sin)
            .mapv_into_any(|v| v < 0.00000001 && v > -0.00000001)
            .fold(true, |a, b| a && *b));
        let true_anomaly_cos = true_anomaly.clone().mapv_into(|v| v.cos());
        let polar_angle_opposite_cos = polar_angle_opposite.clone().mapv_into(|v| v.cos());
        assert!(
            (true_anomaly_cos.mapv_into(|v| v.abs()) - polar_angle_opposite_cos)
                .mapv_into_any(|v| v < 0.00000001 && v > -0.00000001)
                .fold(true, |a, b| a && *b)
        );
    }

    #[test]
    fn azimuthal_angle_for_pi_longitude_of_ascending_node() {
        let argument_of_periapsis: f64 = 0.;
        let inclination: f64 = PI / 4.;
        let true_anomaly: Array1<f64> =
            orbit::calculate_true_anomaly_from_series(Array1::linspace(0., 1., 10), 0., 1., 0.);
        let polar_angle: Array1<f64> = transformations::polar_angle_from_keplerian_elements(
            true_anomaly.clone(),
            inclination,
            argument_of_periapsis,
        );
        let azimuthal_angle_zero_loan: Array1<f64> =
            transformations::azimuthal_angle_from_keplerian_elements_and_polar_angle(
                true_anomaly.clone(),
                polar_angle.clone(),
                argument_of_periapsis,
                0.,
            );
        let azimuthal_angle_pi_loan: Array1<f64> =
            transformations::azimuthal_angle_from_keplerian_elements_and_polar_angle(
                true_anomaly.clone(),
                polar_angle.clone(),
                argument_of_periapsis,
                PI,
            );
        assert_eq!(azimuthal_angle_zero_loan + PI, azimuthal_angle_pi_loan);
    }

    #[test]
    fn azimuthal_angle_for_pi_argument_of_periapsis() {
        let longitude_of_ascending_node: f64 = 0.;
        let inclination: f64 = PI / 4.;
        let argument_of_periapsis: f64 = 0.;
        let half_ticks: usize = 5;
        let ticks: usize = 2 * half_ticks;
        let true_anomaly: Array1<f64> =
            orbit::calculate_true_anomaly_from_series(Array1::linspace(0., 1., ticks), 0., 1., 0.);
        let polar_angle: Array1<f64> = transformations::polar_angle_from_keplerian_elements(
            true_anomaly.clone(),
            inclination,
            argument_of_periapsis,
        );
        let azimuthal_angle_zero_aop: Array1<f64> =
            transformations::azimuthal_angle_from_keplerian_elements_and_polar_angle(
                true_anomaly.clone(),
                polar_angle.clone(),
                0.,
                longitude_of_ascending_node,
            );
        let azimuthal_angle_pi_aop: Array1<f64> =
            transformations::azimuthal_angle_from_keplerian_elements_and_polar_angle(
                true_anomaly.clone(),
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
