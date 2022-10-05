use std::f64::consts::PI;
use std::time::Instant;

fn relax_left(sol: &Vec<f64>, src: &Vec<f64>, s: f64, h: f64) -> f64 {
    sol[1] - src[0] * h.powi(2) / 6.
}

fn relax_middle(sol: &Vec<f64>, src: &Vec<f64>, s: f64, h: f64, i: usize) -> f64 {
    (0. + sol[i + 1] * (1. + h / s) + sol[i - 1] * (1. - h / s) - src[i] * h.powi(2)) / 2.
}

fn relax_right(sol: &Vec<f64>, src: &Vec<f64>, s: f64, h: f64) -> f64 {
    0.
}

fn res_left(sol: &Vec<f64>, src: &Vec<f64>, s: f64, h: f64) -> f64 {
    -(relax_left(sol, src, s, h) - sol[0]) * 6. / h.powi(2)
}

fn res_middle(sol: &Vec<f64>, src: &Vec<f64>, s: f64, h: f64, i: usize) -> f64 {
    -(relax_middle(sol, src, s, h, i) - sol[i]) * 2. / h.powi(2)
}

fn res_right(sol: &Vec<f64>, src: &Vec<f64>, s: f64, h: f64) -> f64 {
    0.
}

const r_s: f64 = 8.;
const rho_c: f64 = 1.28e-3;

fn src(s: &f64) -> f64 {
    if *s < 0.5 {
        let a = 1. - s;
        let b = s / a;
        4. * PI * rho_c * (1. - b.powi(2)) * r_s.powi(2) / a.powi(4)
    } else {
        0.
    }
}

fn sol(s: &f64) -> f64 {
    if *s < 0.5 {
        let a = s / (1. - s);
        -2. * PI * rho_c * r_s.powi(2) * (1. / 2. - a.powi(2) / 3. + a.powi(4) / 10.)
    } else {
        -8. / 15. * PI * rho_c * r_s.powi(2) * (1. - s) / s
    }
}

fn main() {
    let problem = bvp::Problem::new(
        0.,
        1.,
        bvp::Operation::new(relax_left, relax_middle, relax_right),
        bvp::Operation::new(res_left, res_middle, res_right),
        src,
        Some(sol),
    );
    let mut solver = bvp::Solver::new(problem, 16, (4, 1, 4));
    let num_iters = 30;
    let instant = Instant::now();
    for i in 0..num_iters {
        solver.vcycle();
        println!(
            "{} {:e} {:e}",
            i,
            solver.residual_rms_normalized(),
            solver.error_rms()
        );
    }
    println!(
        "averaged time: {}",
        instant.elapsed().as_secs_f64() / (num_iters as f64)
    );
}
