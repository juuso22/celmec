use ndarray::{s, Array1, Array2};

///An orbital mechanics that can be simulated for some time interval and have an impulse applied to it.
///
///Any struct implementing such a system should have at least a starting and ending times for the simulation and a number of time steps to be simulated (sorry, no adaptive time stepping (yet)).
///
///Below, there are instructions on what the methods do when implemented for any struct within `celmec` implementing them. You can use the documentation for inspiration and guidance in case you attempt your own implementations.
pub trait System: Clone {
    ///Get the simulation starting time of a system.
    fn get_start_time(&self) -> f64;
    ///Set the simulation starting time of a system.
    ///
    /// **Inputs**:
    ///
    /// **time**: the new starting time to be set.
    ///
    /// **Output**: A system with the new starting time.
    fn set_start_time(&mut self, time: f64) -> Self;

    ///Get the simulation ending time of a system.
    fn get_end_time(&self) -> f64;
    ///Set the simulation ending time of a system.
    ///
    /// **Inputs**:
    ///
    /// **time**: the new ending time to be set.
    ///
    /// **Output**: A system with the new ending time.
    fn set_end_time(&mut self, time: f64) -> Self;

    ///Get the number of time steps to be simulated for a system.
    fn get_steps(&self) -> usize;
    ///Set the numebr of simulation time steps of a system.
    ///
    /// **Inputs**:
    ///
    /// **steps**: the new number of steps to be set.
    ///
    /// **Output**: A system with the new number of steps.
    fn set_steps(&mut self, steps: usize) -> Self;

    ///Get the simulation step size of a system.
    fn get_step_size(&self) -> f64;
    ///Set the simulation step size from the starting and ending times of the simulation and the number of steps.
    fn set_step_size(&mut self) -> Self;

    ///Set time parameters of the simulation ie. starting time, ending time, number of steps and the step size.
    ///
    /// **Inputs**:
    ///
    /// **start_time**: the new starting time to be set.
    ///
    /// **end_time**: the new ending time to be set.
    ///
    /// **steps**: the new number of steps to be set.
    ///
    /// **Output**: A system with new time parameters.
    fn set_simulation_time(&mut self, start_time: f64, end_time: f64, steps: usize) -> Self;

    ///Get the initial position of a system.
    fn get_initial_rr(&self) -> Array1<f64>;
    ///Get the initial velocity of a system.
    fn get_initial_vv(&self) -> Array1<f64>;
    ///Set the initial conditions (ie. initial position and velocity) of a system.
    ///
    /// **Inputs**:
    ///
    /// **rr0**: the new initial position to be set.
    ///
    /// **vv0**: the new initial velocity to be set.
    ///
    /// **Output**: A system with new initial conditions.
    fn set_initial_conditions(&self, rr0: Array1<f64>, vv0: Array1<f64>) -> Self;

    ///Modify a system by applying an impulse to it.
    ///
    /// **Inputs**:
    ///
    /// **impulse**: an array whose elements are the time of the impulse and its v<sub>x</sub>, v<sub>x</sub> and v<sub>x</sub> components, respectively.
    ///
    /// **output**: a system modified by the impulse
    fn apply_impulse(&self, impulse: Array1<f64>) -> Self;

    ///Do a simulation from the last simulated point to the next by adding the effect of an impulse in-between
    ///
    /// **Inputs**:
    ///
    /// **impulse**: an array whose elements are the time of the impulse and its v<sub>x</sub>, v<sub>x</sub> and v<sub>x</sub> components, respectively.
    ///
    /// **last_simulated_point**: an array whose elements are the time, x, y, z, v<sub>x</sub>, v<sub>x</sub> and v<sub>x</sub> coordinates, respectively.
    ///
    /// **Ouputs**: A system that has its initial conditions set at the next time point after `last_simulated_point` and the impulse applied in between
    ///
    /// To construct the output system, a simulation is first done from `last_simulated_point` to the impulse point using `simulate()`. Then the impulse is applied applied using `apply_impulse()`. If the impulse point does not coincide with a simulation point, `simulate()` is called again to simulate from the impulse point to the next simulation point. Finally, the output system is configured using the final point obtained after all the necessary simulations described above.
    fn simulate_impulse(&self, impulse: Array1<f64>, last_simulated_point: Array1<f64>) -> Self {
        let mut system = self.clone();
        let mut tmp_system = self.clone();

        tmp_system = tmp_system.set_simulation_time(last_simulated_point[0], impulse[0], 2);

        let rr0: Array1<f64> = last_simulated_point.slice(s![1..4]).to_owned();
        let vv0: Array1<f64> = last_simulated_point.slice(s![4..7]).to_owned();
        tmp_system = tmp_system.set_initial_conditions(rr0, vv0);

        let mut patch_orbit: Array2<f64> = self.simulate();
        tmp_system = tmp_system.set_initial_conditions(
            patch_orbit.column(1).slice(s![1..4]).to_owned(),
            patch_orbit.column(1).slice(s![4..7]).to_owned(),
        );
        tmp_system = tmp_system.apply_impulse(impulse.clone());

        let next_start_time: f64 = tmp_system.get_start_time() + self.get_step_size();
        if next_start_time != impulse[0] {
            tmp_system = tmp_system.set_simulation_time(impulse[0], next_start_time, 2);
            patch_orbit = tmp_system.simulate();
            tmp_system = tmp_system.set_initial_conditions(
                patch_orbit.column(1).slice(s![1..4]).to_owned(),
                patch_orbit.column(1).slice(s![4..7]).to_owned(),
            );
        }
        system =
            system.set_initial_conditions(tmp_system.get_initial_rr(), tmp_system.get_initial_vv());
        let new_start_time: f64 = tmp_system.get_end_time();
        let new_end_time: f64 = self.get_end_time();
        system.set_simulation_time(
            new_start_time,
            new_end_time,
            ((new_end_time - new_start_time) / self.get_step_size() + 1.) as usize,
        )
    }

    ///Simulate a system without taking impulses into account.
    ///
    /// **Output**: an array of arrays where is column has the following elements: time, x, y, z, v<sub>x</sub>, v<sub>x</sub> and v<sub>x</sub> coordinates, respectively.
    fn simulate(&self) -> Array2<f64>;
}
