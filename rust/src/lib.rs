mod grid;

#[derive(Clone)]
pub struct Operation {
    left: fn(&Vec<f64>, &Vec<f64>, f64, f64) -> f64,
    middle: fn(&Vec<f64>, &Vec<f64>, f64, f64, usize) -> f64,
    right: fn(&Vec<f64>, &Vec<f64>, f64, f64) -> f64,
}

impl Operation {
    pub fn new(
        left: fn(&Vec<f64>, &Vec<f64>, f64, f64) -> f64,
        middle: fn(&Vec<f64>, &Vec<f64>, f64, f64, usize) -> f64,
        right: fn(&Vec<f64>, &Vec<f64>, f64, f64) -> f64,
    ) -> Self {
        Operation {
            left,
            middle,
            right,
        }
    }
}

#[derive(Clone)]
pub struct Problem {
    range: grid::Range,
    relax: Operation,
    residual: Operation,
    source: fn(&f64) -> f64,
    solution: Option<fn(&f64) -> f64>,
}

impl Problem {
    pub fn new(
        x1: f64,
        x2: f64,
        relax: Operation,
        residual: Operation,
        source: fn(&f64) -> f64,
        solution: Option<fn(&f64) -> f64>,
    ) -> Self {
        Problem {
            range: grid::Range::new(x1, x2),
            relax,
            residual,
            source,
            solution,
        }
    }
}

pub struct Solver {
    problem: Problem,
    solution: grid::Grid,
    source: grid::Grid,
    residual: grid::Grid,
    num_iters: (u32, u32, u32),
}

impl Solver {
    fn new_with_source(
        problem: Problem,
        depth: u32,
        num_iters: (u32, u32, u32),
        source: grid::Grid,
    ) -> Self {
        let range = problem.range.clone();
        Solver {
            problem,
            solution: grid::Grid::new_zeros(range.clone(), depth),
            source,
            residual: grid::Grid::new_zeros(range.clone(), depth),
            num_iters,
        }
    }

    pub fn new(problem: Problem, depth: u32, num_iters: (u32, u32, u32)) -> Self {
        let range = problem.range.clone();
        let source = problem.source;
        Solver::new_with_source(
            problem,
            depth,
            num_iters,
            grid::Grid::new_func(range.clone(), depth, source),
        )
    }
    
    pub fn relax(&mut self) {
        let h = self.solution.coord.h;
        let x = &(self.solution.coord.x);
        let n = self.solution.coord.n;
        let sol_val = &mut (self.solution.val);
        let src_val = &(self.source.val);
        let relax = &(self.problem.relax);
        // red sweep
        sol_val[0] = (relax.left)(sol_val, src_val, x[0], h);
        for i in 1..(n / 2) {
            let j = 2 * i;
            sol_val[j] = (relax.middle)(sol_val, src_val, x[j], h, j);
        }
        sol_val[n] = (relax.right)(sol_val, src_val, x[n], h);
        // black sweep
        for i in 0..(n / 2) {
            let j = 2 * i + 1;
            sol_val[j] = (relax.middle)(sol_val, src_val, x[j], h, j);
        }
    }

    pub fn get_residual(&mut self) {
        let h = self.solution.coord.h;
        let x = &(self.solution.coord.x);
        let n = self.solution.coord.n;
        let sol_val = &(self.solution.val);
        let src_val = &(self.source.val);
        let res_val = &mut (self.residual.val);
        let residual = &(self.problem.residual);
        res_val[0] = (residual.left)(sol_val, src_val, x[0], h);
        for i in 1..n {
            res_val[i] = (residual.middle)(sol_val, src_val, x[i], h, i);
        }
        res_val[n] = (residual.right)(sol_val, src_val, x[n], h);
    }

    pub fn vcycle(&mut self) {
        for _j in 0..self.num_iters.0 {
            self.relax();
        }
        if self.solution.depth > 1 {
            for _j in 0..self.num_iters.1 {
                self.get_residual();
                let mut solver = Solver::new_with_source(
                    self.problem.clone(),
                    self.solution.depth - 1,
                    self.num_iters,
                    self.residual.coarsen(),
                );
                solver.vcycle();
                let error = solver.solution.fine();
                for i in 0..(self.solution.coord.n + 1) {
                    self.solution.val[i] += error.val[i];
                }
            }

            for _j in 0..self.num_iters.2 {
                self.relax();
            }
        }
    }

    fn error(&self) -> grid::Grid {
        let solution = self
            .problem
            .solution
            .expect("Function for solution in Problem should be provided.");
        let mut grid = grid::Grid::new_func(
            self.solution.coord.range.clone(),
            self.solution.depth,
            solution,
        );
        for i in 0..(self.solution.coord.n + 1) {
            grid.val[i] -= self.solution.val[i];
        }
        grid
    }

    pub fn residual_rms(&mut self) -> f64 {
        self.get_residual();
        self.residual.rms()
    }

    pub fn error_rms(&self) -> f64 {
        self.error().rms()
    }
}
