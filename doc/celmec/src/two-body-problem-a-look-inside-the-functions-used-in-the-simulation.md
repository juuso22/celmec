# Two-body Problem: A Look Inside the Functions Used in the Simulation

To see all the physics and maths described in the previous subchapter to actually take place, one can look inside the function `orbit::calculate_f_from_initial_rr_and_vv`. It looks like this:

```
pub fn calculate_f_from_initial_rr_and_vv(
    rr: Array1<f64>,
    vv: Array1<f64>,
    mu: f64,
    start_time: f64,
    end_time: f64,
    steps: usize,
) -> Array1<f64> {
    let ee: Array1<f64> = calculate_ee(rr.clone(), vv.clone(), mu);
    let e: f64 = calculate_e(rr.clone(), vv.clone(), mu);
    let initial_f: f64 = calculate_initial_f_from_initial_conditions(rr.clone(), ee, e);
    let initial_eccentric_anomaly: Array1<f64> =
        calculate_eccentric_anomaly_from_f(array![initial_f], e);
    let a: f64 = calculate_a_from_initial_rr_and_vv(rr, vv, mu);
    let n: f64 = calculate_n(mu, a);
    let tau: f64 = calculate_tau(
        0.,
        calculate_mean_anomaly_from_eccentric_anomaly(initial_eccentric_anomaly, e)[0],
        n,
    );

    let t: Array1<f64> = Array1::linspace(start_time, end_time, steps);
    let eccentric_anomaly: Array1<f64> = calculate_eccentric_anomaly_iteratively(
        t.clone(),
        Array1::zeros(steps),
        0.0001,
        100,
        n,
        e,
        tau,
    );
    calculate_f_from_eccentric_anomaly(eccentric_anomaly, e)
}
```

Let's break that down:

1. `ee` is a vector is calculated with `calculate_ee`. Its length gives the eccentricity `e`, calculated by `calculate_e`.
2. True anomaly `f` at the initial vector position `rr` is calculated in `calculate_initial_f_from_initial_conditions` as it is needed for...
3. ..initial eccentric anomaly calculated in `calculate_eccentric_anomaly_from_f` which in turn is needed for...
4. `tau`, the time at which the planet finds itself at the perihelion, nearest to the sun, calculated in `calculate_tau`.
5. An additional quantity is needed to obtain `tau`. It's called `n` and is calculated in `calculate_n`.
6. Then an array `t` for all the time points in the simulation is initialized.
7. Next is the "big thing" of the simulation when eccentric anomaly is calculated for all simulation time points ie. for the array `t` in `calculate_eccentric_anomaly_iteratively`.
8. Finally we obtain true anomalies `f` from the previously obtained eccentric anomalies in `calculate_f_from_eccentric_anomaly` and return that.

All of the above are part of the `orbit` module of `celmec` and one could drill down further into each of the above to see what they actually do, but to keep it short(ish), I'll just show here that `calculate_eccentric_anomaly_iteratively` is where the Newton-Raphson method gets applied like this:

```
pub fn calculate_eccentric_anomaly_iteratively(
    t: Array1<f64>,
    initial_value: Array1<f64>,
    tolerance: f64,
    max_iterations: usize,
    n: f64,
    e: f64,
    tau: f64,
) -> Array1<f64> {
    let parameters: HashMap<&str, f64> = HashMap::from([("n", n), ("e", e), ("tau", tau)]);
    if e > 1. {
        solve_equation_iteratively(
            &hyperbolic_kepler_eq_iterative_step,
            initial_value,
            t,
            parameters,
            tolerance,
            max_iterations,
        )
    } else if e == 1. {
        solve_equation_iteratively(
            &barker_eq_iterative_step,
            initial_value,
            t,
            parameters,
            tolerance,
            max_iterations,
        )
    } else if (e < 1.) && (e >= 0.) {
        solve_equation_iteratively(
            &kepler_eq_iterative_step,
            initial_value,
            t,
            parameters,
            tolerance,
            max_iterations,
        )
    } else {
        panic!("Eccentricity cannot be negative.");
    }
}
```

I'll not go through that function step by step but instead use it to bring up two points. Firstly, note that there is some equation solving happening for `e` (eccentricity) values of 1 and above. Eccentricity of ellipses is always between 0 and 1 (with 0 included and 1 not). So what does it mean for an orbit to have an eccentricity of 1 or above? It turn out not all orbits are closed ellipses but they can also take the form of a parabola (eccentricity 1) ot a byperbola (eccentricity above 1). These orbits will be the subject of the next subchapter so stay tuned.

But first the second point: in the `calculate_eccentric_anomaly_iteratively` function the last word in its name, `iteratively`, might give a hint that there are other ways to calculate the eccentric anomaly instead of Newton-Raphson. And indeed, one could use a series expansion (if one is familiar with such things) by using `calculate_f_from_series`. Note though, that it currently contains so few terms that it gives something even remotely reliable only for near-circle orbit ie. orbits with eccentricity close to 0.

All the functions described above are public, so they can be used in whatever you might fancy to try with `celmec`. To see the full technical Rust documentation of the library you can either:

1. Git clone [the `celmec` repo](https://github.com/juuso22/celmec.git) and run `cargo doc` inside it if you are familiar with git. The output of `cargo doc` should tell you where to look next.
2. Navigate to URL_HERE to see the same docs.

Further in this book, we'll tear some other used functions open in case the author has deemed them to have eaten things of interest. This is done to help out with navigating the the technical Rust documentation as I find those sometimes intimidating to tacḱle head on and want to ease the burden for anyone else feeling similarly.
