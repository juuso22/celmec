use crate::orbital_elements::KeplerianElements;
use ndarray::Array1;
use std::collections::HashMap;

/// A system struct contains information about what is to be simulated: what physics model is to be used, the necessary parameters to carry out simulations using the chosen model
#[derive(Clone)]
pub struct System {
    pub model: String,
    pub start_time: f64,
    pub end_time: f64,
    pub steps: usize,
    pub array_parameters: HashMap<String, Array1<f64>>,
    pub f64_parameters: HashMap<String, f64>,
    pub perturbations: Vec<String>,
}

/// Create a 2-body system
///
/// **Inputs**:
///
/// start_time: time at which the simulation of the system starts
///
/// end_time: time at which the simulation of the system ends
///
/// steps: number of time steps to be simulated
///
/// mu: gravitational parameter of the system
///
/// rr0: position of the system at start_time
///
/// vv0: velocity of the system at start_time
///
/// **Output**: a system struct for 2-body system
pub fn create_2_body_system(
    start_time: f64,
    end_time: f64,
    steps: usize,
    mu: f64,
    rr0: Array1<f64>,
    vv0: Array1<f64>,
) -> System {
    let mut system = System {
        model: "2-body".to_string(),
        start_time,
        end_time,
        steps,
        array_parameters: HashMap::new(),
        f64_parameters: HashMap::new(),
        perturbations: vec![],
    };
    system.f64_parameters.insert("mu".to_string(), mu);
    system.array_parameters.insert("rr0".to_string(), rr0);
    system.array_parameters.insert("vv0".to_string(), vv0);
    system
}

///// Create a 2-body system from keplerian elements
/////
///// **Inputs**:
/////
///// start_time: time at which the simulation of the system starts
/////
///// end_time: time at which the simulation of the system ends
/////
///// steps: number of time steps to be simulated
/////
///// mu: gravitational parameter of the system
/////
///// elements: keplerian elements of the 2-body system
/////
///// **Output**: a system struct for 2-body system
//pub fn create_2_body_system(
//    start_time: f64,
//    end_time: f64,
//    steps: usize,
//    mu: f64,
//    elements: KeplerianElements,
//) -> System {
//    //TODO
//}

/// Validates as system struct.
///
/// **Inputs**:
///
/// system: a system struct
///
/// Specific validation steps are described in the subfunctions.
pub fn validate_system(system: System) {
    if system.model == "2-body".to_string() {
        validate_2_body_system(system);
    }
}

/// Validates a system struct for a 2-body simulation
///
/// **Inputs**:
///
/// system: a system struct
///
/// A 2-body system requires a gravitational parameter mu, and some initial position rr0 and initial velocity vv0.
///
/// If instead of rr0 and vv0 you have keplerian elements for the system. Use [create_2_body_system](`create_2_body_system`) to create a system which has the initial conditions calculated from the elements.
pub fn validate_2_body_system(system: System) {
    assert!(system.f64_parameters.contains_key("mu"));
    assert!(system.array_parameters.contains_key("rr0"));
    assert!(system.array_parameters.contains_key("vv0"));
}
