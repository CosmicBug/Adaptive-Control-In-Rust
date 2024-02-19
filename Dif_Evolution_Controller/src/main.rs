extern crate nalgebra as na;
extern crate csv;
extern crate serde;
extern crate differential_evolution;

extern crate rayon;
use rayon::prelude::*;

use differential_evolution::self_adaptive_de;
use na::Vector4; // For a 4-dimensional vector
use std::f32::consts::PI;
use serde::Serialize;
use std::error::Error;

#[derive(Serialize)]
struct StateRecord {
    t: f32,
    theta: f32,
    phi: f32,
    p_theta: f32,
    p_phi: f32,
}


// Functiom to calculate the derivatives of the state vector - θ, φ, p_θ, p_φ
// i.e equations of motion

fn derivatives(state: &Vector4<f32>, t: f32, v_scalar: f32) -> Vector4<f32> {
    let m = 1.0; // mass
    let l: f32 = 1.0; // length
    let g = 9.81; // acceleration due to gravity

    let (θ, φ, pθ, pφ) = (state[0], state[1], state[2], state[3]);
    let cosθφ = (θ - φ).cos();
    let sinθφ = (θ - φ).sin();
    let cosθ = θ.cos();
    let sinθ = θ.sin();
    let cosφ = φ.cos();
    let sinφ = φ.sin();

    let denominator = m * l.powi(2) * (cosθφ.powi(2) - 2.0);

    let θ_dot = (pφ * cosθφ - pθ - m * l * v_scalar * (cosφ * cosθφ - 2.0 * cosθ)) / denominator;
    let φ_dot = (pθ * cosθφ - 2.0 * pφ - 2.0 * m * l * v_scalar * (cosθ * cosθφ - cosφ)) / denominator;
    let pθ_dot = -2.0 * m * l * v_scalar * θ_dot * sinθ - m * l.powi(2) * θ_dot * φ_dot * sinθφ - 2.0 * m * g * l * sinθ;
    let pφ_dot = -m * l * v_scalar * φ_dot * sinφ + m * l.powi(2) * θ_dot * φ_dot * sinθφ - m * g * l * sinφ;

    Vector4::new(θ_dot, φ_dot, pθ_dot, pφ_dot)
}


// Function to perform a single step of the Runge-Kutta 4th order method
fn rk4_step(f: &dyn Fn(&Vector4<f32>, f32, f32) -> Vector4<f32>, state: &Vector4<f32>, t: f32, dt: f32, v_scalar: f32) -> Vector4<f32> {
    let k1 = dt * f(state, t, v_scalar);
    let k2 = dt * f(&(state + 0.5 * k1), t + 0.5 * dt, v_scalar);
    let k3 = dt * f(&(state + 0.5 * k2), t + 0.5 * dt, v_scalar);
    let k4 = dt * f(&(state + k3), t + dt, v_scalar);

    state + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
}


#[derive(Serialize)]
struct EnergyRecord {
    t: f32,
    e_total: f32,
}

// Function to calculate the total energy of the system
fn calculate_total_energy(m: f32, l: f32, v: f32, θ: f32, φ: f32, θ_dot: f32, φ_dot: f32, g: f32) -> f32 {
    let cosθ = θ.cos();
    let cosφ = φ.cos();
    let cosθφ = (θ - φ).cos();

    0.5 * m * (2.0 * v.powi(2) + 2.0 * l.powi(2) * θ_dot.powi(2) + 4.0 * v * l * θ_dot * cosθ + l.powi(2) * φ_dot.powi(2) + 2.0 * v * l * φ_dot * cosφ + 2.0 * l.powi(2) * θ_dot * φ_dot * cosθφ) - m * g * l * (2.0 * cosθ + cosφ - 3.0)
}




fn main() -> Result<(), Box<dyn Error>> {
    let dt = 0.001;
    let t_final = 25.0;
    let n_steps = (t_final / dt) as usize;

    // Parameters for the INPUT cosine signal
    let a = 1.0; // Amplitude
    let ω = 2.0 * PI / t_final; // Angular frequency
    let φ_signal = 0.0; // Phase shift

    let mut state = Vector4::new(1.5, -0.5, 0.0, 0.0);
    let mut v_scalar = vec![0.0; n_steps];

    let mut wtr = csv::Writer::from_path("state_evolution.csv")?;

    // Simulate the system for n_steps
    for i in 0..n_steps {
        let t = i as f32 * dt;
        v_scalar[i] = a * (ω * t + φ_signal).cos();
        state = rk4_step(&derivatives, &state, t, dt, v_scalar[i]);

        // Create a record for the current state
        let record = StateRecord {
            t,
            theta: state[0],
            phi: state[1],
            p_theta: state[2],
            p_phi: state[3],
        };

        // Write the record to the CSV
        // wtr.serialize(record)?;
    }

    let m = 1.0; // mass
    let l = 1.0; // length
    let g = 9.81; // acceleration due to gravity

    let mut energy_writer = csv::Writer::from_path("energy_evolution.csv")?;

    for i in 0..n_steps {
        let t = i as f32 * dt;
        v_scalar[i] = a * (ω * t + φ_signal).cos(); // Update v_scalar for each timestep

        let current_v = v_scalar[i];
        let state = rk4_step(&derivatives, &state, t, dt, current_v);
        let θ = state[0];
        let φ = state[1];
        let pθ = state[2];
        let pφ = state[3];
        // Derive θ_dot and φ_dot from the state, assuming derivatives function returns them
        let derivatives_result = derivatives(&state, t, current_v);
        let θ_dot = derivatives_result[0];
        let φ_dot = derivatives_result[1];

        let e_total = calculate_total_energy(m, l, current_v, θ, φ, θ_dot, φ_dot, g);

        // energy_writer.serialize(EnergyRecord { t, e_total })?;

        if i == n_steps - 1 {
            println!("t = {:.3}, E_total = {:.3}", t, e_total);
        }
    }


    let mut wtr_opt = csv::Writer::from_path("optimized_control_input.csv")?;

    let search_area = vec![(-0.5, 0.5); n_steps];

    let mut de = self_adaptive_de(search_area, move |v: &[f32]| {
        // let mut state = Vector4::new(1.5, -0.5, 0.0, 0.0);
        let initial_state = Vector4::new(0.2, 0.3, 0.0, 0.0); // Initial conditions of the system
        let m = 1.0;
        let l = 1.0;
        let g = 9.81;
        let mut e_total_sum = 0.0;

        // for i in 0..n_steps {
        // Parallel computation of the total energy sum
        let e_total_sum: f32 = (0..n_steps).into_par_iter().map(|i| {
            let t = i as f32 * dt;
            let current_v = v[i];
            let mut state = initial_state;
            state = rk4_step(&derivatives, &state, t, dt, current_v);
            let derivatives_result = derivatives(&state, t, current_v);
            // let e_total = calculate_total_energy(m, l, g, current_v, &state, &derivatives_result);
            let θ = state[0];
            let φ = state[1];
            let pθ = state[2];
            let pφ = state[3];
            // Derive θ_dot and φ_dot from the state, assuming derivatives function returns them
            let derivatives_result = derivatives(&state, t, current_v);
            let θ_dot = derivatives_result[0];
            let φ_dot = derivatives_result[1];

            // let e_total = calculate_total_energy(m, l, current_v, θ, φ, θ_dot, φ_dot, g);

            // e_total_sum += e_total;
            calculate_total_energy(m, l, current_v, θ, φ, θ_dot, φ_dot, g)
        }).sum();

        e_total_sum
        // e_total_sum
    });

    de.iter().nth(10000);

    let (cost, pos) = de.best().unwrap();
    println!("Minimum energy sum: {}", cost);


        // Write the optimized v vector to the CSV file
    // Assuming pos is a slice &[f32]
    for &value in pos.iter() {
        wtr_opt.serialize(value)?;
    }

    // Ensure all data is flushed to the file
    wtr_opt.flush()?;
    // println!("Optimized v vector: {:?}", pos);


    wtr.flush()?;
    energy_writer.flush()?;
    Ok(())
}