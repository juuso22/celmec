/// Systems module contains a system struct which represents what is to be simulated and validation functions for that struct to ensure for all necessary information is included for a given system.
pub mod system;

///This module contains "well known" constants
pub mod constants;

/// This module contains mathematical helper functions for orbit calculations.
pub mod math;

/// This module contains functions that deal with mechanics more generally than within the confines of the celestial variant
pub mod mechanics;

/// This module contains functions make calculations of about orbits of 2-body systems
///
/// ## Naming
/// Naming tends to cause headaches at times and this module is no different. Here are some terms I use in a way that are to my best knowledge non-standard, but for which the standard terminology escapes me:
///
/// Eccentric anomaly: I call the equivalent of the eccentric anomaly for elliptic orbit eccentric anomaly also in the context of prabolic and hyperbolic orbits.
///
/// 'hyperbolic Kepler's equation': This is the equation M = E - e * sinh(E) which relates mean anomaly M to eccentric anomaly (see note above) for hyperbolic orbits (with help from eccentricity e). I don't know what the correct name for this equation is.
pub mod two_body;

/// This module contains structs for orbital elements.
///
/// It is currently limited to deal with the 2-body problem.
pub mod orbital_elements;

///This module contains functions to apply impulses to orbit
pub mod impulse;

/// This module contains functions for transformations between coordinate systems.
///
/// Eg. switching from the orbital plane of a 2-body system to spherical coordinates with one of the bodies at origin.
pub mod transformations;

#[cfg(test)]
mod math_tests {
    use super::*;
    use math::*;
    use ndarray::{array, Array1};
    use std::collections::HashMap;
    use std::f64::consts::PI;

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

    #[test]
    fn atan() {
        assert!(math::atan(0., 0.).is_nan());
        assert_eq!(math::atan(1., 0.), PI / 2.);
        assert_eq!(math::atan(-1., 0.), -PI / 2.);
        assert_eq!(math::atan(1. / 2_f64.sqrt(), 1. / 2_f64.sqrt()), PI / 4.);
        assert_eq!(
            math::atan(1. / 2_f64.sqrt(), -1. / 2_f64.sqrt()),
            3. * PI / 4.
        );
        assert_eq!(math::atan(-1. / 2_f64.sqrt(), 1. / 2_f64.sqrt()), -PI / 4.);
        assert_eq!(
            math::atan(-1. / 2_f64.sqrt(), -1. / 2_f64.sqrt()),
            5. * PI / 4.
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{array, Array1};

    #[test]
    fn k_dot_e_is_almost_zero() {
        let rr0: Array1<f64> = array![1., 3.2, 0.4];
        let v0: Array1<f64> = array![-4., 0.002, 2.2];
        let mu: f64 = 5.123;
        let ee: Array1<f64> = orbital_elements::calculate_ee(rr0.clone(), v0.clone(), mu);
        let kk: Array1<f64> = math::cross_product(rr0, v0);
        assert!((ee * kk).sum() < 1.0e-10)
    }
}

#[cfg(test)]
mod two_body_tests {
    use super::*;
    use ndarray::{array, Array1, Array2};
    use orbital_elements::*;
    use std::f64::consts::PI;
    use two_body::*;

    #[test]
    fn simple_mean_anomaly() {
        assert_eq!(
            array![-4., -5.],
            calculate_mean_anomaly(array![1., 0.5], 2., 3.)
        );
    }

    #[test]
    fn simple_f_from_series() {
        assert_eq!(
            array![PI, PI / 2. + 1.75 - 13. / 12.],
            calculate_f_from_series(array![1., 0.5], 1., 2., 0.)
        );
    }

    #[test]
    fn r_from_f_with_zero_e() {
        let a: f64 = 1.;
        assert_eq!(array![a], calculate_r_from_f(array![3.], 0., a));
    }

    #[test]
    fn velocity() {
        let rr0: Array1<f64> = array![1., 0., 0.];
        let mu: f64 = 2.;

        //Elliptic case
        let vv0: Array1<f64> = array![0., 1. / 2_f64.sqrt(), 1. / 2_f64.sqrt()];
        let e: f64 = calculate_e(rr0.clone(), vv0.clone(), mu);
        let a: f64 = calculate_a_from_initial_rr_and_vv(rr0.clone(), vv0.clone(), mu);
        let eccentric_anomaly: Array1<f64> = calculate_eccentric_anomaly_from_initial_rr_and_vv(
            rr0.clone(),
            vv0.clone(),
            mu,
            0.,
            0.,
            1,
        );
        let v: Array1<f64> = calculate_v_from_eccentric_anomaly(eccentric_anomaly, mu, e, a);
        assert!((math::euclidean_norm(vv0.clone()) - v[0]) < 1.0e-10);
        let vv: Array2<f64> = calculate_vv_from_v_rr_and_initial_vv(
            v,
            array![[rr0[0]], [rr0[1]], [rr0[2]]],
            vv0.clone(),
        );
        assert!((vv0[0] - vv.column(0)[0]) < 1.0e-10);
        assert!((vv0[1] - vv.column(0)[1]) < 1.0e-10);
        assert!((vv0[2] - vv.column(0)[2]) < 1.0e-10);

        //Hyperbolic case
        let vv0: Array1<f64> = array![0., 3. / 2_f64.sqrt(), 1. / 2_f64.sqrt()];
        let e: f64 = calculate_e(rr0.clone(), vv0.clone(), mu);
        let a: f64 = calculate_a_from_initial_rr_and_vv(rr0.clone(), vv0.clone(), mu);
        let eccentric_anomaly: Array1<f64> = calculate_eccentric_anomaly_from_initial_rr_and_vv(
            rr0.clone(),
            vv0.clone(),
            mu,
            0.,
            0.,
            1,
        );
        let v: Array1<f64> = calculate_v_from_eccentric_anomaly(eccentric_anomaly, mu, e, a);
        assert!((math::euclidean_norm(vv0.clone()) - v[0]) < 1.0e-10);
        let vv: Array2<f64> = calculate_vv_from_v_rr_and_initial_vv(
            v,
            array![[rr0[0]], [rr0[1]], [rr0[2]]],
            vv0.clone(),
        );
        assert!((vv0[0] - vv.column(0)[0]) < 1.0e-10);
        assert!((vv0[1] - vv.column(0)[1]) < 1.0e-10);
        assert!((vv0[2] - vv.column(0)[2]) < 1.0e-10);
    }
}

#[cfg(test)]
mod transformations_tests {
    use super::*;
    use ndarray::{array, Array1, Array2};
    use orbital_elements::KeplerianElements;
    use rand::prelude::*;
    use std::f64::consts::PI;
    use transformations::*;

    #[test]
    fn phi_for_zero_iota() {
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
                - phi_from_keplerian_elements(
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
                - phi_from_keplerian_elements(
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
        let phi: Array1<f64> =
            phi_from_keplerian_elements(f.clone(), iota, omega, longitude_of_the_ascending_node);
        assert_eq!(
            ((f.clone() + omega).mapv_into(|v| v.sin()) - phi.clone().mapv_into(|v| v.sin()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            f.clone().len() as u8
        );
        assert_eq!(
            ((f.clone() + omega).mapv_into(|v| v.cos()) - phi.mapv_into(|v| v.cos()))
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
                - phi_from_keplerian_elements(
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
    fn theta_for_zero_iota() {
        let iota: f64 = 0.;
        let mut omega: f64 = 0.;
        let f: Array1<f64> = array![0., PI / 2., PI, 3. / 2. * PI];
        assert_eq!(
            array![PI / 2., PI / 2., PI / 2., PI / 2.],
            theta_from_keplerian_elements(f.clone(), iota, omega,)
        );

        let mut rng = rand::thread_rng();
        omega = rng.gen::<f64>() * 2. * PI;
        println!("Random omega: {}", omega);
        assert_eq!(
            array![PI / 2., PI / 2., PI / 2., PI / 2.],
            theta_from_keplerian_elements(f, iota, omega,)
        );
    }

    #[test]
    fn phi_for_pi_per_two_iota() {
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
        let phi: Array1<f64> =
            phi_from_keplerian_elements(f.clone(), iota, omega, longitude_of_the_ascending_node);
        assert_eq!(
            (array![0., 0., 0., 0., 0., 0., 0., 0.] - phi.clone().mapv_into(|v| v.sin()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            6
        );
        assert_eq!(
            (array![1., 1., -1., -1., -1., 1., 1., 1.] - phi.clone().mapv_into(|v| v.cos()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            6
        );
        assert_eq!(phi.mapv_into_any(|v| (v.is_nan()) as u8).sum(), 2);

        let mut rng = rand::thread_rng();
        longitude_of_the_ascending_node = rng.gen::<f64>() * 2. * PI;
        println!(
            "Random longitude of the ascending node: {}",
            longitude_of_the_ascending_node
        );
        let phi: Array1<f64> =
            phi_from_keplerian_elements(f.clone(), iota, omega, longitude_of_the_ascending_node);
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
        assert_eq!(
            (control.clone().mapv_into(|v| v.sin()) - phi.clone().mapv_into(|v| v.sin()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            6
        );
        assert_eq!(
            (control.clone().mapv_into(|v| v.cos()) - phi.clone().mapv_into(|v| v.cos()))
                .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
                .sum(),
            6
        );
        assert_eq!(phi.mapv_into_any(|v| (v.is_nan()) as u8).sum(), 2);
    }

    #[test]
    fn theta_for_pi_per_two_iota() {
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
                    v - PI / 2.
                } else if (v >= PI) & (v < PI * 3. / 2.) {
                    v - PI / 2.
                } else if v >= PI * 3. / 2. {
                    5. * PI / 2. - v
                } else {
                    PI / 2. - v
                }
            }) - theta_from_keplerian_elements(f.clone(), iota, omega,))
            .mapv_into_any(|v| (v.abs() < 1.0e-10) as u8)
            .sum(),
            f.clone().len() as u8
        );

        f = array![PI];
        assert!(0. - theta_from_keplerian_elements(f, iota, omega,)[0] < 1.0e-10);
    }

    #[test]
    fn cartesians_from_sphericals() {
        let theta: Array1<f64> = array![PI / 4.];
        let phi: Array1<f64> = array![PI / 2.];
        let r: Array1<f64> = array![2_f64.sqrt()];
        let xyz: Array2<f64> = cartesian_coordinates_from_spherical_coordinates(phi, theta, r);
        assert!((xyz[[0, 0]] - 0.) < 1.0e-10);
        assert!((xyz[[1, 0]] - 1.) < 1.0e-10);
        assert!((xyz[[2, 0]] - 1.) < 1.0e-10);
    }

    #[test]
    fn cartesians_from_polars() {
        let theta: Array1<f64> = array![3. * PI / 4.];
        let r: Array1<f64> = array![2. * 2_f64.sqrt()];
        let xyz: Array2<f64> = cartesian_coordinates_from_polar_coordinates(theta, r);
        assert!((xyz[[0, 0]] + 2.) < 1.0e-10);
        assert!((xyz[[1, 0]] - 2.) < 1.0e-10);
        assert_eq!(xyz[[2, 0]], 0.);
    }

    #[test]
    fn sphericals_for_halleys_comet() {
        //This is more like a regression test:
        //I picked the values by adapting the demo coordinate_transformation
        //at a time when I was pretty sure the functions tested here work correctly
        //(ie. 16.10.2024)
        //The modification I made was to calculate spherical coordinates instead of
        //cartesians.
        let iota: f64 = 2.82674;
        let omega: f64 = 1.95564;
        let longitude_of_the_ascending_node: f64 = 1.03666;
        let f: Array1<f64> = array![
            3.1396486577792015,
            -3.0618489166515466,
            -2.907447721707956,
            2.8766713188022917,
            3.0513809874254925,
        ];
        assert_eq!(
            array![
                1.8621521944438066,
                1.8513135484528696,
                1.8257708767582101,
                1.8833107655695367,
                1.8716463184907006,
            ],
            theta_from_keplerian_elements(f.clone(), iota, omega,)
        );
        assert_eq!(
            array![
                6.149735192706617,
                6.234129679841356,
                0.10886093708854638,
                5.875107345304375,
                6.057983613318866,
            ],
            phi_from_keplerian_elements(f.clone(), iota, omega, longitude_of_the_ascending_node)
        );
    }

    #[test]
    fn sphericals_from_f_and_r() {
        let keplerian_elements = KeplerianElements {
            e: 0.,
            a: 1.,
            tau: 0.,
            iota: 0.0,
            omega: PI / 2.,
            longitude_of_the_ascending_node: f64::NAN,
        };
        let mu: f64 = 1.;
        let f: Array1<f64> =
            two_body::calculate_f_from_keplerian_elements(&keplerian_elements, mu, 0., 2. * PI, 9);
        println!("True anomaly f: {:#?}", f.clone());

        let sphericals: Array2<f64> =
            spherical_coordinates_from_f_and_keplerian_elements(f.clone(), keplerian_elements);
        //Theta tests
        assert_eq!(sphericals[[1, 0]], PI / 2.);
        assert!(
            (sphericals.row(1).to_owned() - Array1::<f64>::ones(9) * PI / 2.)
                .iter()
                .sum::<f64>()
                .abs()
                < 1.0e-10
        );
        //Phi tests
        println!("Phi: {:#?}", sphericals.row(0));
        assert!(
            (sphericals.row(0).to_owned()
                - array![
                    2. * PI / 4.,
                    3. * PI / 4.,
                    4. * PI / 4.,
                    5. * PI / 4.,
                    6. * PI / 4.,
                    7. * PI / 4.,
                    8. * PI / 4.,
                    1. * PI / 4.,
                    2. * PI / 4.,
                ])
            .iter()
            .sum::<f64>()
            .abs()
                < 1.0e-10
        );
    }

    #[test]
    fn cartesians_from_f_and_r() {
        let keplerian_elements = KeplerianElements {
            e: 0.,
            a: 1.,
            tau: 0.,
            iota: 0.0,
            omega: PI / 2.,
            longitude_of_the_ascending_node: f64::NAN,
        };
        let mu: f64 = 1.;
        let f: Array1<f64> =
            two_body::calculate_f_from_keplerian_elements(&keplerian_elements, mu, 0., 3., 30);
        println!("True anomaly f: {:#?}", f.clone());
        let r: Array1<f64> =
            two_body::calculate_r_from_f(f.clone(), keplerian_elements.e, keplerian_elements.a);
        println!("Radius r: {:#?}", r.clone());

        let xyz: Array2<f64> = cartesian_coordinates_from_f_r_and_keplerian_elements(
            f.clone(),
            r.clone(),
            keplerian_elements,
        );
        assert_eq!(xyz.row(2).into_iter().sum::<f64>(), 0.);
        assert!(
            (xyz.row(0).map(|v| v.powf(2.)) + xyz.row(1).map(|v| v.powf(2.))
                - r.mapv_into(|v| v.powf(2.)))
            .iter()
            .sum::<f64>()
            .abs()
                < 1.0e-10
        );
    }
}

#[cfg(test)]
mod orbital_elements_tests {
    use super::*;
    use math::*;
    use ndarray::{array, Array1};
    use orbital_elements::*;
    use std::f64::consts::PI;
    use two_body::*;

    #[test]
    fn iota_from_initial_conditions() {
        let rr0: Array1<f64> = array![1., 0., 0.];
        let vv0: Array1<f64> = array![0., 1., 0.];
        assert_eq!(calculate_iota_from_initial_rr_and_vv(rr0, vv0), 0.);
        let rr0: Array1<f64> = array![1., 0., 0.];
        let vv0: Array1<f64> = array![0., 0., 1.];
        assert_eq!(calculate_iota_from_initial_rr_and_vv(rr0, vv0), PI / 2.);
        let rr0: Array1<f64> = array![0., 1. / 2_f64.sqrt(), 1. / 2_f64.sqrt()];
        let vv0: Array1<f64> = array![-1., 0., 0.];
        assert_eq!(calculate_iota_from_initial_rr_and_vv(rr0, vv0), PI / 4.);
        let rr0: Array1<f64> = array![0., -1. / 2_f64.sqrt(), 1. / 2_f64.sqrt()];
        let vv0: Array1<f64> = array![-1., 0., 0.];
        assert_eq!(
            calculate_iota_from_initial_rr_and_vv(rr0, vv0),
            3. * PI / 4.
        );
        let rr0: Array1<f64> = array![0., -1. / 2_f64.sqrt(), -1. / 2_f64.sqrt()];
        let vv0: Array1<f64> = array![1., 0., 0.];
        assert_eq!(calculate_iota_from_initial_rr_and_vv(rr0, vv0), PI / 4.);
        let rr0: Array1<f64> = array![0., 1. / 2_f64.sqrt(), -1. / 2_f64.sqrt()];
        let vv0: Array1<f64> = array![1., 0., 0.];
        assert_eq!(
            calculate_iota_from_initial_rr_and_vv(rr0, vv0),
            3. * PI / 4.
        );
    }

    #[test]
    fn longitude_of_the_ascending_node_from_initial_conditions() {
        let rr0: Array1<f64> = array![1., 0., 0.];
        let vv0: Array1<f64> = array![0., 1., 0.];
        assert!(
            calculate_longitude_of_the_ascending_node_from_initial_rr_and_vv(rr0, vv0).is_nan()
        );
        let rr0: Array1<f64> = array![1., 0., 0.];
        let vv0: Array1<f64> = array![0., 1., 1.];
        assert_eq!(
            calculate_longitude_of_the_ascending_node_from_initial_rr_and_vv(rr0, vv0),
            0.
        );
        let rr0: Array1<f64> = array![0., 1., 0.];
        let vv0: Array1<f64> = array![-1., 0., 1.];
        assert_eq!(
            calculate_longitude_of_the_ascending_node_from_initial_rr_and_vv(rr0, vv0),
            PI / 2.
        );
        let rr0: Array1<f64> = array![0., 1., 0.];
        let vv0: Array1<f64> = array![1., 0., 1.];
        assert_eq!(
            calculate_longitude_of_the_ascending_node_from_initial_rr_and_vv(rr0, vv0),
            PI / 2.
        );
        //TODO: more cases
    }

    #[test]
    fn earth_e() {
        let solar_mass: f64 = 1.98847 * 10_f64.powf(30.);
        let earth_mass: f64 = 5.972 * 10_f64.powf(24.);
        let mu: f64 = two_body::calculate_mu(solar_mass, earth_mass);
        let earth_perihelion = array![147.10 * 10_f64.powf(9.), 0., 0.];
        let earth_orbital_vv_at_perihelion = array![0., 30290., 0.];
        assert!(
            calculate_e(
                earth_perihelion.clone(),
                earth_orbital_vv_at_perihelion.clone(),
                mu
            ) - 0.0167
                < 0.0005
        );
        assert!(
            calculate_e(earth_perihelion, earth_orbital_vv_at_perihelion, mu) - 0.0167 > -0.0005
        );
    }

    #[test]
    fn a_and_e_for_circular_orbit() {
        let rr0: Array1<f64> = array![0., -1., 0.];
        let vv0: Array1<f64> = array![0., 0., 1.];
        let mu: f64 = 1.;
        let e: f64 = calculate_e(rr0.clone(), vv0.clone(), mu);
        assert_eq!(e, 0.);
        let h: f64 = calculate_h(rr0.clone(), vv0, mu);
        let a: f64 = calculate_a(mu, h);
        assert_eq!(a, euclidean_norm(rr0));
    }

    #[test]
    fn omega() {
        let rr0: Array1<f64> = array![0., -1., 0.];
        let vv0: Array1<f64> = array![0., 0., 1.];
        let mu: f64 = 1.;
        //Circular orbit, so cannot determine perihelion
        assert!(calculate_omega_from_mu_and_initial_rr_and_vv(mu, rr0, vv0).is_nan());

        //Venus: https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
        let rr0: Array1<f64> = array![
            107.480e6 * 0.08204_f64.cos(),
            0.,
            107.480e6 * 0.08204_f64.sin()
        ];
        let vv0: Array1<f64> = array![0., -35.26e3, 0.];
        let mu: f64 = calculate_mu(1.989e30, 4.8673e24);
        assert_eq!(
            calculate_omega_from_mu_and_initial_rr_and_vv(mu, rr0, vv0),
            PI / 2.
        );

        //TODO: come up with more cases
    }
}
