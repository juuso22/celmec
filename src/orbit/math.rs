use ndarray::{array, Array1};
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
