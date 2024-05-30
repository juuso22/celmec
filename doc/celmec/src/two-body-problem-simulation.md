# Two Body Problem: Simulation

As a first simulation with `celmec`, we'll simulate one round of Mercury around the sun as a two body problem. We'll first walk through building the simulation, but the final Rust code can be found further down.

## Walking through the writing of the simulation

First init a new Rust project called `two-body-problem`:

```
cargo init two-body-problem
```

Edit the project's `Cargo.toml` by adding the following under `[dependencies]`:

```
ndarray = "0.15.6" #Or choose your preferred version
celmec = { git = "https://github.com/juuso22/celmec.git" }
```

Then we can start editing the project's `src/main.rs`. First, remove the default content and import the following from `ndarray` 

```
use ndarray::{Array, Array1}
```

Then import the following from celmec:

```
use celmec::orbit;
```

The `orbit` module contains functions functions for the simulation.

Next, we need the masses of Mercury and the sun and some initial conditions for the position and velocity of the planet with respect to the sun. Sun's mass is (the googlable) 1.989 * 10<sup>30</sup>kg and the mass of Mercury (along other info about the planet that will be used), namely 3.3010 * 10<sup>23</sup> kg, can be found from [Nasa's fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html). To have some simple initial conditions, let's have Mercury furthest away from the sun it reaches (ie. Mercury's aphelion) and let's fix our coordinates so that this point is to the direction of the positive x-axis. The distance of Mercury from the sun in the aphelion is 69.818 * 10<sup>9</sup> m. According to Kepler's third law, Mercury's orbital velocity is the slowest at aphelion. That slowest speed is 38.86 * 10<sup>3</sup> m/s and that is in the direction of y-axis. Let's decide, for convenience, that is the to the direction of the positive y-axis. Then we write this all into the `main` function of our project:

```
    let mu: f64 = orbit::calculate_mu(3.301e23, 1.989e30);
    let r0: Array1<f64> = array![69.818e9, 0., 0.];
    let v0: Array1<f64> = array![0., 38860., 0.];
```



## Final Rust code

`Cargo.toml`:

```
[package]
name = "two-body-problem"
version = "0.1.0"
edition = "2021"

[dependencies]
ndarray = "0.15.6"
celmec = { git = "https://github.com/juuso22/celmec.git" }
```

`main.rs`:

```
use celmec::orbit;
use ndarray::{array, Array, Array1};
use std::fs::File;
use std::io::Write;

fn main() {
    let mu: f64 = orbit::calculate_mu(3.301e23, 1.989e30);
    let r0: Array1<f64> = array![69.818e9, 0., 0.];
    let v0: Array1<f64> = array![0., -38860., 0.];

    let ee: Array1<f64> = orbit::calculate_ee(r0.clone(), v0.clone(), mu);
    let e: f64 = orbit::calculate_e(r0.clone(), v0.clone(), mu);
    let h: f64 = orbit::calculate_h(r0.clone(), v0, mu);
    let a: f64 = orbit::calculate_a(mu, h);
    let n: f64 = orbit::calculate_n(mu, a);
    let initial_f: f64 = orbit::calculate_initial_f_from_initial_conditions(r0, ee, e);
    let initial_eccentric_anomaly: Array1<f64> =
        orbit::calculate_eccentric_anomaly_from_f(array![initial_f], e);
    let tau: f64 = orbit::calculate_tau(
        0.,
        orbit::calculate_mean_anomaly_from_eccentric_anomaly(initial_eccentric_anomaly, e)[0],
        n,
    );

    let ticks: usize = 100;
    let t: Array1<f64> = Array::linspace(0., 7603200., ticks);

    let eccentric_anomaly: Array1<f64> = orbit::calculate_eccentric_anomaly_iteratively(
        t.clone(),
        Array::zeros(ticks),
        0.0001,
        100,
        orbit::calculate_n(mu, a),
        e,
        tau,
    );
    let f: Array1<f64> = orbit::calculate_f_from_eccentric_anomaly(eccentric_anomaly, e);
    let radius: Array1<f64> = orbit::calculate_radius_from_f(f.clone(), e, a);

    println!("Eccentricity: {}", e);
    println!(
        "Maximum distance from the sun (aphelion): {}",
        radius.iter().max_by(|a, b| a.total_cmp(b)).unwrap()
    );
    println!(
        "Minimum distance from the sun (perihelion): {}",
        radius.iter().min_by(|a, b| a.total_cmp(b)).unwrap()
    );
    println!("Perihelion time: {}", tau / 3600. / 24.);

    let mut coordinate_file = File::create("mercury.csv").unwrap();
    write!(coordinate_file, "t,radius,f\n").unwrap();
    for i in 0..=(ticks - 1) {
        write!(coordinate_file, "{},{},{}\n", t[i], radius[i], f[i]).unwrap();
    }
}
```
