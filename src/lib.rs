/// This module contains functions to calculate various properties of orbits.
///
/// It is currently limited to deal with the 2-body problem.
///
/// ## Naming
/// Naming tends to cause headaches at times and this module is no different. Here are some terms I use in a way that are to my best knowledge non-standard, but for which the standard terminology escapes me:
///
/// Eccentric anomaly: I call the equivalent of the eccentric anomaly for elliptic orbit eccentric anomaly also in the context of prabolic and hyperbolic orbits.
///
/// 'hyperbolic Kepler's equation': This is the equation M = E - e * sinh(E) which relates mean anomaly M to eccentric anomaly (see note above) for hyperbolic orbits (with help from eccentricity e). I don't know what the correct name for this equation is.
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
    use rand::prelude::*;
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
    fn r_from_f_with_zero_e() {
        let a: f64 = 1.;
        assert_eq!(array![a], orbit::calculate_r_from_f(array![3.], 0., a));
    }

    #[test]
    fn theta_for_zero_iota() {
        let iota: f64 = 0.;
        let mut omega: f64 = 0.;
        let mut longitude_of_the_ascending_node: f64 = 0.;
        let f: Array1<f64> = array![
            PI / 4.,
            PI / 2.,
            PI * 3. / 4.,
            PI,
            PI * 5. / 4.,
            PI * 3. / 2.,
            PI * 7. / 4.,
            0.0
        ];

        assert_eq!(
            (f.clone()
                - transformations::theta_from_keplerian_elements(
                    f.clone(),
                    iota,
                    omega,
                    longitude_of_the_ascending_node
                ))
            .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
            .sum(),
            f.clone().len() as u8
        );

        longitude_of_the_ascending_node = PI / 4.;
        assert_eq!(
            (f.clone() + longitude_of_the_ascending_node
                - transformations::theta_from_keplerian_elements(
                    f.clone(),
                    iota,
                    omega,
                    longitude_of_the_ascending_node
                ))
            .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
            .sum(),
            f.clone().len() as u8
        );

        longitude_of_the_ascending_node = 0.;
        omega = PI / 4.;
        //Here a funny corner case happens:
        //Due to numerical inaccuracies, there is a case where the values to be compared are 0.0
        //and something a tiny bit less than 2*PI. Therefore, the previously used check of a small
        //difference in the angles themselves does not work out and we compare sines and cosines
        //of the relevant angles instead.
        let theta: Array1<f64> = transformations::theta_from_keplerian_elements(
            f.clone(),
            iota,
            omega,
            longitude_of_the_ascending_node,
        );
        assert_eq!(
            ((f.clone() + omega).mapv_into(|v| v.sin()) - theta.clone().mapv_into(|v| v.sin()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            f.clone().len() as u8
        );
        assert_eq!(
            ((f.clone() + omega).mapv_into(|v| v.cos()) - theta.mapv_into(|v| v.cos()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            f.clone().len() as u8
        );

        let mut rng = rand::thread_rng();
        longitude_of_the_ascending_node = rng.gen::<f64>() * 2. * PI;
        println!(
            "Random longitude of the ascending node: {}",
            longitude_of_the_ascending_node
        );
        omega = rng.gen::<f64>() * 2. * PI;
        println!("Random argument of perihelion (omega): {}", omega);
        assert_eq!(
            ((f.clone() + omega + longitude_of_the_ascending_node) % (2. * PI)
                - transformations::theta_from_keplerian_elements(
                    f.clone(),
                    iota,
                    omega,
                    longitude_of_the_ascending_node
                ))
            .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
            .sum(),
            f.clone().len() as u8
        );
    }

    #[test]
    fn phi_for_zero_iota() {
        let iota: f64 = 0.;
        let mut omega: f64 = 0.;
        let f: Array1<f64> = array![0., PI / 2., PI, 3. / 2. * PI];
        assert_eq!(
            array![0., 0., 0., 0.],
            transformations::phi_from_keplerian_elements(f.clone(), iota, omega,)
        );

        let mut rng = rand::thread_rng();
        omega = rng.gen::<f64>() * 2. * PI;
        println!("Random omega: {}", omega);
        assert_eq!(
            array![0., 0., 0., 0.],
            transformations::phi_from_keplerian_elements(f, iota, omega,)
        );
    }

    #[test]
    fn theta_for_pi_per_two_iota() {
        let iota: f64 = PI / 2.;
        let omega: f64 = 0.;
        let mut longitude_of_the_ascending_node: f64 = 0.;
        let f: Array1<f64> = array![
            PI / 4.,
            PI / 2.,
            PI * 3. / 4.,
            PI,
            PI * 5. / 4.,
            PI * 3. / 2.,
            PI * 7. / 4.,
            0.0
        ];

        //Here a funny corner case happens:
        //Due to numerical inaccuracies, there is a case where the values to be compared are 0.0
        //and something a tiny bit less than 2*PI. Therefore, the previously used check of a small
        //difference in the angles themselves does not work out and we compare sines and cosines
        //of the relevant angles instead.
        let theta: Array1<f64> = transformations::theta_from_keplerian_elements(
            f.clone(),
            iota,
            omega,
            longitude_of_the_ascending_node,
        );
        assert_eq!(
            (array![0., 0., 0., 0., 0., 0., 0., 0.] - theta.clone().mapv_into(|v| v.sin()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            f.clone().len() as u8
        );
        assert_eq!(
            (array![1., 1., -1., -1., -1., 1., 1., 1.] - theta.mapv_into(|v| v.cos()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            f.clone().len() as u8
        );

        let mut rng = rand::thread_rng();
        longitude_of_the_ascending_node = rng.gen::<f64>() * 2. * PI;
        println!(
            "Random longitude of the ascending node: {}",
            longitude_of_the_ascending_node
        );
        let theta: Array1<f64> = transformations::theta_from_keplerian_elements(
            f.clone(),
            iota,
            omega,
            longitude_of_the_ascending_node,
        );
        let control: Array1<f64> = array![
            longitude_of_the_ascending_node,
            longitude_of_the_ascending_node,
            (PI + longitude_of_the_ascending_node) % (2. * PI),
            (PI + longitude_of_the_ascending_node) % (2. * PI),
            (PI + longitude_of_the_ascending_node) % (2. * PI),
            longitude_of_the_ascending_node,
            longitude_of_the_ascending_node,
            longitude_of_the_ascending_node
        ];
        println!("{}", theta);
        println!("{}", control);
        assert_eq!(
            (control.clone().mapv_into(|v| v.sin()) - theta.clone().mapv_into(|v| v.sin()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            f.clone().len() as u8
        );
        assert_eq!(
            (control.clone().mapv_into(|v| v.cos()) - theta.mapv_into(|v| v.cos()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            f.clone().len() as u8
        );
    }

    #[test]
    fn phi_for_pi_per_two_iota() {
        let iota: f64 = PI / 2.;
        let omega: f64 = 0.;
        let mut f: Array1<f64> = array![
            PI / 5.,
            PI * 2. / 5.,
            PI * 3. / 5.,
            PI * 4. / 5.,
            PI * 6. / 5.,
            PI * 7. / 5.,
            PI * 8. / 5.,
            PI * 9. / 5.,
            0.0
        ];
        assert_eq!(
            (f.clone().mapv_into(|v| {
                if (v >= PI / 2.) & (v < PI) {
                    PI - v
                } else if (v >= PI) & (v < PI * 3. / 2.) {
                    -v + PI
                } else if v >= PI * 3. / 2. {
                    v - 2. * PI
                } else {
                    v
                }
            }) - transformations::phi_from_keplerian_elements(f.clone(), iota, omega,))
            .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
            .sum(),
            f.clone().len() as u8
        );

        f = array![PI];
        assert!(0. - transformations::phi_from_keplerian_elements(f, iota, omega,)[0] < 1.0e-10);
    }
}
