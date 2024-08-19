# Keplerian Elements: Coordinate Transformation

So far all the visualisations in this book have been restricted to the orbital plane of the system being simulated. However, the Keplerian elements contain a bunch of angles, which determine the the position of the orbit with respect to the orbital plane of the Earth. For 3d visuals, let us know use these angles to switch from the plane of the system to a three-dimensional volume where the Sun lies at origin and the xy-plane is the orbital plane of the earth and use spherical coordinates there. In order to do that, we can copy the code from the original simulation of Halley's comet's orbit and add some pieces.

Firstly, in order to do the coordinate switch, we'll need the `transformations` module of `celmec` so we'll modify our imports as follows:

```
use celmec::{orbit, orbital_elements, transformations};
```

Then, after calculating the radii and true anomalies, the new module can be used to perform the coordinate transformation:

```
    let spherical_coordinates: (Array1<f64>, Array1<f64>) =
        transformations::spherical_coordinates_from_f_and_keplerian_elements(
            f,
            halleys_keplerian_elements,
        );

```

Then the output file generation has to be adjusted as well. Instead of what we previously had, we'll use:

```
    let mut coordinate_file = File::create("halleys_orbit_3d.csv").unwrap();
    write!(coordinate_file, "t,radius,theta,phi\n").unwrap();
    for i in 0..=(ticks - 1) {
        write!(
            coordinate_file,
            "{},{},{},{}\n",
            time[i], radius[i], spherical_coordinates.0[i], spherical_coordinates.1[i]
        )
        .unwrap();
    }
```

### The Final Rust Code

In the end, `Config.toml` should look like this:

```
[package]
name = "coordinate-transformation"
version = "0.1.0"
edition = "2021"

[dependencies]
ndarray = "0.15.6"
celmec = { git = "https://github.com/juuso22/celmec.git" }
```

and `main.rs` should look like this

```
use celmec::{orbit, orbital_elements, transformations};
use ndarray::{Array, Array1};
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;

fn main() {
    let halleys_keplerian_elements = orbital_elements::KeplerianElements {
        e: 0.96658,
        longitude_of_the_ascending_node: 1.03666,
        tau: 0.,
        a: 2.65342e12,
        iota: 2.82674,
        omega: 1.95564,
    };

    let ticks = 50;
    let rotation_time: f64 = 2379801600.;
    let time: Array1<f64> = Array::linspace(
        halleys_keplerian_elements.tau - rotation_time / 2.,
        halleys_keplerian_elements.tau + rotation_time / 2.,
        ticks,
    );
    let eccentric_anomaly: Array1<f64> = orbit::calculate_eccentric_anomaly_iteratively(
        time.clone(),
        time.clone(),
        0.00001,
        100,
        2. * PI / rotation_time,
        halleys_keplerian_elements.e,
        halleys_keplerian_elements.tau,
    );

    let f: Array1<f64> =
        orbit::calculate_f_from_eccentric_anomaly(eccentric_anomaly, halleys_keplerian_elements.e);

    let radius = orbit::calculate_r_from_f(
        f.clone(),
        halleys_keplerian_elements.e,
        halleys_keplerian_elements.a,
    );

    let spherical_coordinates: (Array1<f64>, Array1<f64>) =
        transformations::spherical_coordinates_from_f_and_keplerian_elements(
            f,
            halleys_keplerian_elements,
        );

    let mut coordinate_file = File::create("halleys_orbit_3d.csv").unwrap();
    write!(coordinate_file, "t,radius,theta,phi\n").unwrap();
    for i in 0..=(ticks - 1) {
        write!(
            coordinate_file,
            "{},{},{},{}\n",
            time[i], radius[i], spherical_coordinates.0[i], spherical_coordinates.1[i]
        )
        .unwrap();
    }
}
```

### Visualising the results

For the visualisation, the following python script can be used:

```
TODO: add script
```

which can be run like this:

```
TODO
```

This should give something like this:

TODO: pic and gif here
