use ndarray::{s, Array1, Array2};

///TODO: what should go into trait docs?
pub trait System: Clone {
    //fn clone(&self) -> Self;

    fn get_start_time(&self) -> f64;
    fn set_start_time(&mut self, time: f64) -> Self;

    fn get_end_time(&self) -> f64;
    fn set_end_time(&mut self, time: f64) -> Self;

    fn get_steps(&self) -> usize;
    fn set_steps(&mut self, steps: usize) -> Self;

    fn get_step_size(&self) -> f64;
    fn set_step_size(&mut self) -> Self;

    fn set_simulation_time(&mut self, start_time: f64, end_time: f64, steps: usize) -> Self;

    //TODO: make more general
    fn set_initial_conditions(&self, rr0: Array1<f64>, vv0: Array1<f64>) -> Self;
    fn get_initial_rr(&self) -> Array1<f64>;
    fn get_initial_vv(&self) -> Array1<f64>;

    fn apply_impulse(&self, impulse: Array1<f64>) -> Self;

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

    fn simulate(&self) -> Array2<f64>;
}
