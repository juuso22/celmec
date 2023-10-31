use ndarray::{array, Array1};
use std::collections::HashMap;
//use std::f64::sqrt;

pub fn cross_product(a: Array1<f64>, b: Array1<f64>) -> Array1<f64> {
    let c1: f64 = a[1] * b[2] - a[2] * b[1];
    let c2: f64 = -a[0] * b[2] + a[2] * b[0];
    let c3: f64 = a[0] * b[1] - a[1] * b[0];
    array![c1, c2, c3]
}

pub fn euclidean_norm(a: Array1<f64>) -> f64 {
    (a[0].powf(2.) + a[1].powf(2.) + a[2].powf(2.)).sqrt()
}

pub fn solve_equation_iteratively(
    f: &dyn Fn(Array1<f64>, Array1<f64>, HashMap<&str, f64>) -> Array1<f64>,
    initial_value: Array1<f64>,
    time: Array1<f64>,
    parameters: HashMap<&str, f64>,
    tolerance: f64,
    max_iterations: usize,
) -> Array1<f64> {
    let mut old_value: Array1<f64> = initial_value;
    let mut new_value: Array1<f64> = f(old_value.clone(), time.clone(), parameters.clone());
    let mut iteration_round: usize = 1;
    while ((old_value.clone() - new_value.clone())
        .mapv_into(|v| v.abs())
        .fold(0_f64, |m, v| m.max(*v))
        > tolerance)
        && (iteration_round < max_iterations)
    {
        old_value = new_value.clone();
        new_value = f(old_value.clone(), time.clone(), parameters.clone());
        iteration_round += 1;
    }
    new_value
}
