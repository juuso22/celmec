use crate::math::cross_product;
use ndarray::Array1;

/// Calculates angular momentum per unit of mass for a 2 body problem from initial conditions.
///
/// **Inputs**:
///
/// rr: The initial position of the rotating body with respect to the central body.
///
/// vv: The initial velocity of the rotating body with respect to the central body.
///
/// **Output**: angular momentum per unit of mass.
pub fn calculate_kk_from_initial_rr_and_vv(rr: Array1<f64>, vv: Array1<f64>) -> Array1<f64> {
    cross_product(rr, vv)
}
