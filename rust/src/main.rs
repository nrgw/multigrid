use std::f64::consts::PI;
use std::time::Instant;

fn relax_left(sol: &[f64], src: &[f64], s: f64, h: f64) -> f64 {
    sol[1] - src[0] * h.powi(2) / 6.
}

fn relax_middle(sol: &[f64], src: &[f64], s: f64, h: f64, i: usize) -> f64 {
    (0. + sol[i + 1] * (1. + h / s) + sol[i - 1] * (1. - h / s) - src[i] * h.powi(2)) / 2.
}

fn relax_right(sol: &[f64], src: &[f64], s: f64, h: f64) -> f64 {
    0.
}

fn res_left(sol: &[f64], src: &[f64], s: f64, h: f64) -> f64 {
    -(relax_left(sol, src, s, h) - sol[0]) * 6. / h.powi(2)
}

fn res_middle(sol: &[f64], src: &[f64], s: f64, h: f64, i: usize) -> f64 {
    -(relax_middle(sol, src, s, h, i) - sol[i]) * 2. / h.powi(2)
}

fn res_right(sol: &[f64], src: &[f64], s: f64, h: f64) -> f64 {
    0.
}

const K: f64 = 100.;
const N: f64 = 1.;
const rho_c: f64 = 1.28e-3;

trait Sinc {
    fn sinc(&self) -> f64;
}

impl Sinc for f64 {
    fn sinc(&self) -> f64 {
        if *self == 0. {
            1.
        } else {
            self.sin() / self
        }
    }
}

fn src(s: &f64) -> f64 {
    if *s < 0.5 {
        let a = ((N + 1.) * K * rho_c.powf(1. / N  - 1.) / (4. * PI)).sqrt();
        let r_s = PI * a;    
        let rho = rho_c * (PI * s / (1. - s)).sinc();
        4. * PI * rho * r_s.powi(2) / (1. - s).powi(4)
    } else {
        0.
    }
}

fn sol(s: &f64) -> f64 {
    let a_sq = (N + 1.) * K * rho_c.powf(1. / N  - 1.) / (4. * PI);
    let factor = if *s < 0.5 {
        1. + (PI * s / (1. - s)).sinc()
    } else {
        (1. - s) / s
    };

    -4. * PI * rho_c * a_sq * factor      
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
