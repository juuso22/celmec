# Keplerian Elements: A Look Inside the Transformation

To wrap up with Keplerian elements, let's have a look inside the coordinate transformation function `cartesian_coordinates_from_f_r_and_keplerian_elements`:

```
pub fn cartesian_coordinates_from_f_r_and_keplerian_elements(
    f: Array1<f64>,
    r: Array1<f64>,
    keplerian_elements: orbital_elements::KeplerianElements,
) -> (Array1<f64>, Array1<f64>, Array1<f64>) {
    let (theta, phi) = spherical_coordinates_from_f_and_keplerian_elements(f, keplerian_elements);
    cartesian_coordinates_from_spherical_coordinates(theta, phi, r)
}
```

There we see that to obtain the cartesinan coordinates, we actually go through spherical coordinates first. Let's dig deeper to see where those come from in `spherical_coordinates_from_f_and_keplerian_elements`:

```
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
```

Here are separate functions for both spherical angle coordinates `theta` and `phi`. Both functions to obtain these cooedinates, `theta_from_keplerian_elements` and `phi_from_keplerian_elements`, are public for you to use as you please if using `celmec`. The technical Rust documentation will in turn provide details about the formulae used to calculate the angles: go check it out if you want to tear open the black box of maths used by `celmec`!
