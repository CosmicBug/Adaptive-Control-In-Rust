
extern crate nalgebra as na;
#[allow(unused_imports)]
use na::{SMatrix, DMatrix, Matrix, LU}; // 
#[allow(unused_imports)]
use na::{SVector, DVector, Vector, Vector6}; //
#[allow(unused_imports)]
use na::{U6, U1};
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;


// states: x, x_dot, theta, p_theta, phi, p_phi

struct LQRController {
    q: SMatrix<f64, 6, 6>,
    r: SMatrix<f64, 1, 1>,
    k: SMatrix<f64, 1, 6>,
}

impl LQRController {
    fn new() -> Result<LQRController, &'static str> {
        let q = SMatrix::<f64, 6, 6>::identity();
        let r = SMatrix::<f64, 1, 1>::identity();
        let k = SMatrix::<f64, 1, 6>::zeros();

        Ok(LQRController {
            q,
            r,
            k,
        })
    }

    #[allow(dead_code)] // this is available only for discrete-time scheme, i.e., x_(k+1) = A * x_k + B * u_k
    fn calculate_gain_dare(&mut self,
        a: &SMatrix<f64, 6, 6>,
        b: &SMatrix<f64, 6, 1>,
        q: &SMatrix<f64, 6, 6>,
        r: &SMatrix<f64, 1, 1>,
        tol: f64
    ) -> Result<SMatrix<f64, 1, 6>, &'static str> { // credits from Ryohei Sasaki, https://github.com/rsasaki0109/rust_robotics/blob/main/src/bin/lqr_steer_control.rs
        let mut x = q.clone(); 
        let max_iter = 1000;
        // let tol = 0.01;
        let mut iter = 0;
        while iter < max_iter {
            let xn = a.transpose() * x * a - a.transpose() * x * b *
                (r + b.transpose() * x * b).try_inverse().expect("Couldn't compute inverse")
                * b.transpose() * x * a + q;
    
            if (xn - x).abs().max() < tol {
                break;
            }
            x = xn;
            iter += 1;
        }

        self.k = (b.transpose() * x * b + r).try_inverse().expect("Couldn't compute inverse") * (b.transpose() * x * a);
        // self.k = k;
        Ok(self.k.clone())

    }

    fn _p_eom(&self, _t: f64, 
        p: &SMatrix<f64, 6, 6>, 
        a: &SMatrix<f64, 6, 6>,
        b: &SMatrix<f64, 6, 1>,
        q: &SMatrix<f64, 6, 6>,
        r_inv: &SMatrix<f64, 1, 1>,
    ) -> SMatrix<f64, 6, 6> {
        -(a.transpose() * p + p * a - p * b * r_inv * b.transpose() * p + q)
    }
    fn _p_rk4(&self, t: f64, h: f64,
        p: &SMatrix<f64, 6, 6>,
        a: &SMatrix<f64, 6, 6>,
        b: &SMatrix<f64, 6, 1>,
        q: &SMatrix<f64, 6, 6>,
        r_inv: &SMatrix<f64, 1, 1>
    ) -> SMatrix<f64, 6, 6> {
        let k1 = self._p_eom(t, p, a, b, q, r_inv);
        let k2 = self._p_eom(t + h / 2.0, &(p + h / 2.0 * k1), a, b, q, r_inv);
        let k3 = self._p_eom(t + h / 2.0, &(p + h / 2.0 * k2), a, b, q, r_inv);
        let k4 = self._p_eom(t + h, &(p + h * k3), a, b, q, r_inv);
        p + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    }

    fn calculate_gain_care(&mut self,
        a: &SMatrix<f64, 6, 6>,
        b: &SMatrix<f64, 6, 1>,
        q: &SMatrix<f64, 6, 6>,
        r: &SMatrix<f64, 1, 1>,
        tol: f64,
        t_step: f64
    ) -> Result<SMatrix<f64, 1, 6>, &'static str> {
        // calculate the gain matrix using the continuous-time algebraic Riccati equation (CARE) and Runge-Kutta 4th order method with reverse time integration
        // let t_step = t_step_init; // set initial time step (0.01)
        let mut t = 0.0;
        let r_inv = r.try_inverse().ok_or("Couldn't compute inverse of R")?;
        let mut p = q.clone(); // set initial p matrix
        let max_iter = 100000;
        let mut iter = 0;
        while iter < max_iter {
            // let dp_dt = -(a.transpose() * p + p * a - p * b * &r_inv * b.transpose() * p + q);

            // let p_prev = p - t_step * dp_dt; // update p according to the time step (Euler's method, not stable enough)
            let p_prev = self._p_rk4(t, -t_step, &p, a, b, q, &r_inv);
            
            if (p_prev - p).norm() < tol {
                break;
            }

            p = p_prev;
            iter += 1;
            t -= t_step;
        }

        if iter == max_iter {
            println!("Max iterations reached, result may not be accurate!");
        }

        self.k = &r_inv * b.transpose() * p;
    
        let k = self.k.clone();
        let sys_control_mat = a - b * k;
        let eigenvals = sys_control_mat.schur().complex_eigenvalues();
        let eigenvals_real: Vec<f64> = eigenvals.iter().map(|c| c.re).collect();
        println!("Eigenvalues:\n{:?}", eigenvals_real);
        for val in eigenvals_real {
            if val > 0.0 {
                return Err("Unstable system!");
            }
        }
    
        Ok(k)

    }
}

struct DoublePendulum {
    m: f64,                // mass of each pendulum
    g: f64,                // acceleration due to gravity
    l: f64,                // length of each pendulum
    s: f64,                // sign of the pendulum (1 for normal, -1 for inverted)
    a: SMatrix<f64, 6, 6>, // original system matrix (Ax)
    b: SMatrix<f64, 6, 1>, // control input matrix   (Bu)
    q: SMatrix<f64, 6, 6>, // LQR cost matrix (x^T Q x)
    r: SMatrix<f64, 1, 1>, // LQR cost matrix (u^T R u)
    k: SMatrix<f64, 1, 6>, // LQR gain matrix
    c: SMatrix<f64, 6, 6>, // control system matrix (A - BK)
}

impl DoublePendulum {
    fn new(m: f64, g: f64, l: f64, s: f64) -> DoublePendulum {
        let a = SMatrix::<f64, 6, 6>::from_row_slice(&[
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, -1.0 / l * s, 0.0, 1.0 / (m * l * l), 0.0, -1.0 / (m * l * l),
            0.0, 0.0, -2.0 * m * g * l * s, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -1.0 / (m * l * l), 0.0, 2.0 / (m * l * l),
            0.0, 0.0, 0.0, 0.0, -m * g * l * s, 0.0,
        ]);
        let b = SMatrix::<f64, 6, 1>::from_column_slice(&[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);

        DoublePendulum {
            m, g, l, s, a, b,
            q: SMatrix::<f64, 6, 6>::identity(),
            r: SMatrix::<f64, 1, 1>::identity(),
            k: SMatrix::<f64, 1, 6>::zeros(),
            c: SMatrix::<f64, 6, 6>::zeros(),
        }
    }

    fn check_controlability(&self) -> bool {
        let a = self.a.clone();
        let mut b = self.b.clone();
        let n = self.a.nrows();
        let k = self.b.ncols();
        let mut c = DMatrix::<f64>::zeros(n, n*k);
        for j in 0..k {
            c.set_column(j, &b.column(j));
        }
        for i in 1..n {
            b = a * b;
            // println!("B_col: {:?}", b);
            for j in 0..k {
                c.set_column(i*k+j, &b.column(j));
            }
        }
        // println!("Ctrb Matrix:\n{:?}", c);
        let u = c.lu().u();
        let mut rank_c = 0;
        for i in 0..n {
            if u[(i, i)].abs() > 1e-10 {
                rank_c += 1;
            }
        }
        assert!(rank_c <= n, "Rank of ctrb_mat should not > n");

        rank_c == n

        // c.determinant() != 0.0
    }

    fn add_lqr(&mut self, q: &SMatrix<f64, 6, 6>, r: &SMatrix<f64, 1, 1>) {
        self.q = *q;
        self.r = *r;
    }

    fn add_gain(&mut self) {
        let mut controller = LQRController::new().unwrap();
        let result = controller.calculate_gain_care(&self.a, &self.b, &self.q, &self.r, 1e-6, 0.01);
        let k = match result {
            Ok(k) => k,
            Err(e) => panic!("Error: {}", e),
        };
        self.k = k;
        self.c = self.a.clone() - self.b.clone() * self.k.clone();
    }

    fn equations_of_motion(&self, _t: f64, state: Vector6<f64>, active: bool) -> Vector6<f64> {
        let (_x, x_dot, theta, p_theta, phi, p_phi) = (
            state[0], state[1], state[2], state[3], state[4], state[5],
        );

        let v = x_dot;
        let mut acc = 0.0;
        if active {
            if self.s == 1.0 {
                acc = (self.c * state)[1];
            }
            else {
                acc = (self.c * (state - Vector6::new(0., 0., PI, 0., PI, 0.)))[1];
            }
        }

        let theta_dot = (p_phi * (theta - phi).cos() - p_theta - self.m * self.l * v * ((phi).cos() * (theta - phi).cos() - 2.0 * (theta).cos()))
            / ((self.m * self.l.powi(2)) * ((theta - phi).cos().powi(2) - 2.0));
        let phi_dot = (p_theta * (theta - phi).cos() - 2.0 * p_phi - 2.0 * self.m * self.l * v * ((theta).cos() * (theta - phi).cos() - (phi).cos()))
            / ((self.m * self.l.powi(2)) * ((theta - phi).cos().powi(2) - 2.0));

        let p_theta_dot = -2.0 * self.m * self.l * v * theta_dot * (theta).sin()
            - self.m * self.l.powi(2) * theta_dot * phi_dot * (theta - phi).sin()
            - 2.0 * self.m * self.g * self.l * (theta).sin();
        let p_phi_dot = -self.m * self.l * v * phi_dot * (phi).sin()
            + self.m * self.l.powi(2) * theta_dot * phi_dot * (theta - phi).sin()
            - self.m * self.g * self.l * (phi).sin();

        Vector6::new(v, acc, theta_dot, p_theta_dot, phi_dot, p_phi_dot)
    }

    fn rk4(&self, t: f64, state: Vector6<f64>, h: f64, active: bool) -> Vector6<f64> {
        let k1 = self.equations_of_motion(t, state, active);
        let k2 = self.equations_of_motion(t + h / 2.0, state + h / 2.0 * k1, active);
        let k3 = self.equations_of_motion(t + h / 2.0, state + h / 2.0 * k2, active);
        let k4 = self.equations_of_motion(t + h, state + h * k3, active);
        state + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    }

    fn tot_energy(&self, state: Vector6<f64>) -> f64 {
        let (_x, x_dot, theta, p_theta, phi, p_phi) = (
            state[0], state[1], state[2], state[3], state[4], state[5],
        );
        let theta_dot = (p_phi * (theta - phi).cos() - p_theta - self.m * self.l * x_dot * ((phi).cos() * (theta - phi).cos() - 2.0 * (theta).cos()))
            / ((self.m * self.l.powi(2)) * ((theta - phi).cos().powi(2) - 2.0));
        let phi_dot = (p_theta * (theta - phi).cos() - 2.0 * p_phi - 2.0 * self.m * self.l * x_dot * ((theta).cos() * (theta - phi).cos() - (phi).cos()))
            / ((self.m * self.l.powi(2)) * ((theta - phi).cos().powi(2) - 2.0));
        
        0.5 * self.m * (
            2.0 * x_dot * x_dot +
            2.0 * self.l * self.l * theta_dot * theta_dot +
            4.0 * x_dot * self.l * theta_dot * theta.cos() +
            self.l * self.l * phi_dot * phi_dot +
            2.0 * x_dot * self.l * phi_dot * phi.cos() +
            2.0 * self.l * self.l * theta_dot * phi_dot * (theta - phi).cos()
        ) - self.m * self.g * self.l * (
            2.0 * theta.cos() + phi.cos() - 3.0
        )
    }

    fn _k_eom(&self, _t: f64, state: Vector6<f64>, k: &SMatrix<f64, 1, 6>, active: bool) -> Vector6<f64> {
        (self.a.clone() - self.b.clone() * k) * state
    }
    fn _k_rk4(&self, t: f64, h: f64, state: Vector6<f64>, k: &SMatrix<f64, 1, 6>, active: bool) -> Vector6<f64> {
        let k1 = self._k_eom(t, state, k, active);
        let k2 = self._k_eom(t + h / 2.0, state + h / 2.0 * k1, k, active);
        let k3 = self._k_eom(t + h / 2.0, state + h / 2.0 * k2, k, active);
        let k4 = self._k_eom(t + h, state + h * k3, k, active);
        state + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    }
    fn delta_hamiltonian(&self,
        x: f64, x_dot: f64, theta: f64, p_theta: f64, phi: f64, p_phi: f64,
        k: &SMatrix<f64, 1, 6>, t_step: f64
    ) -> f64 {
        let state = Vector6::new(x, x_dot, theta, p_theta, phi, p_phi);
        let t = 0.0; // set current time to 0
        let state_next = self._k_rk4(t, t_step, state, k, true);
        self.tot_energy(state_next) - self.tot_energy(state)
    }

    fn simulate(&self, initial_state: Vector6<f64>, t_end: f64, t_step: f64, active: bool) -> (Vec<Vector6<f64>>, Vec<f64>, Vec<f64>) {
        // prepare the state and energy vectors
        let mut states = vec![];
        let mut energies = vec![];
        let mut times = vec![];
        // initial state
        let mut t = 0.0;
        let mut state = initial_state;
        states.push(state);
        energies.push(self.tot_energy(state));
        times.push(t);
        // simulate the system
        while t < t_end {
            state = self.rk4(t, state, t_step, active);
            states.push(state);
            energies.push(self.tot_energy(state));
            times.push(t);
            t += t_step;
        }
        (states, energies, times)
    }
}


fn main() {
    // define the system variables
    let m:f64 = 1.0;
    let g:f64 = 10.0;
    let l:f64 = 1.0;
    // define the LQR cost matrices
    // 1., 1., 1., 1., 1., 1.
    // 1., 10., 50., 100., 50., 100.
    // 10.
    let q_vals = Vector6::new(1., 10., 50., 100., 50., 100.);
    let q = SMatrix::<f64, 6, 6>::from_diagonal(&q_vals);
    let r = SMatrix::<f64, 1, 1>::from_diagonal_element(10.);
    // define t_end and t_step
    let t_end = 20.0;
    let t_step = 0.01;

    // create the double pendulum system
    let mut double_pendulum = DoublePendulum::new(m, g, l, 1.0);
    let initial_state = Vector6::new(0., 0., 0.2, 0., 0.3, 0.);
    // check the controlability of the system
    println!("Controlability: {}", double_pendulum.check_controlability());

    // simulate the double pendulum system without control
    let (states, energies, times) = double_pendulum.simulate(initial_state, t_end, t_step, false);

    // write the results to a file
    let mut file = File::create("double_pendulum_no_control.csv").unwrap();
    writeln!(file, "t,x,x_dot,theta,p_theta,phi,p_phi,energy").unwrap();
    for (i, state) in states.iter().enumerate() {
        writeln!(file, "{},{},{},{},{},{},{},{}", times[i], state[0], state[1], state[2], state[3], state[4], state[5], energies[i]).unwrap();
    }

    // add lqr into the system
    double_pendulum.add_lqr(&q, &r);
    double_pendulum.add_gain();
    println!("LQR Gain: {:?}", double_pendulum.k);
    // simulate the system with control
    let (states, energies, times) = double_pendulum.simulate(initial_state, t_end, t_step, true);

    // write the results to a file
    let mut file = File::create("double_pendulum.csv").unwrap();
    writeln!(file, "t,x,x_dot,theta,p_theta,phi,p_phi,energy").unwrap();
    for (i, state) in states.iter().enumerate() {
        writeln!(file, "{},{},{},{},{},{},{},{}", times[i], state[0], state[1], state[2], state[3], state[4], state[5], energies[i]).unwrap();
    }

    println!("------------------------------------------");

    // construct the inverted pendulum system
    // note that the balance point is at theta=PI, phi=PI
    let mut inverted_double_pendulum = DoublePendulum::new(m, g, l, -1.0);
    let initial_state = Vector6::new(0., 0., PI - 0.2, 0., PI - 0.3, 0.);
    println!("Controlability (inverted): {}", inverted_double_pendulum.check_controlability());

    // simulate the inverted double pendulum system without control
    let (states, energies, times) = inverted_double_pendulum.simulate(initial_state, t_end, t_step, false);

    // write the results to a file
    let mut file = File::create("inverted_double_pendulum_no_control.csv").unwrap();
    writeln!(file, "t,x,x_dot,theta,p_theta,phi,p_phi,energy").unwrap();
    for (i, state) in states.iter().enumerate() {
        writeln!(file, "{},{},{},{},{},{},{},{}", times[i], state[0], state[1], state[2], state[3], state[4], state[5], energies[i]).unwrap();
    }

    // add lqr into the system
    inverted_double_pendulum.add_lqr(&q, &r);
    inverted_double_pendulum.add_gain();
    println!("LQR Gain (inverted): {:?}", inverted_double_pendulum.k);
    // simulate the system
    let (states, energies, times) = inverted_double_pendulum.simulate(initial_state, t_end, t_step, true);

    // write the results to a file
    let mut file = File::create("inverted_double_pendulum.csv").unwrap();
    writeln!(file, "t,x,x_dot,theta,p_theta,phi,p_phi,energy").unwrap();
    for (i, state) in states.iter().enumerate() {
        writeln!(file, "{},{},{},{},{},{},{},{}", times[i], state[0], state[1], state[2], state[3], state[4], state[5], energies[i]).unwrap();
    }
    
}
