#[derive(Debug, Clone)]
pub struct Range(f64, f64);

impl Range {
    pub fn new(x1: f64, x2: f64) -> Self {
        assert!(x1 < x2, "x1: {} should be smaller than x2: {}.", x1, x2);
        Range(x1, x2)
    }

    fn tuple(&self) -> (f64, f64) {
        (self.0, self.1)
    }
}

#[derive(Debug)]
pub struct Coord {
    pub range: Range,
    pub n: usize,
    pub h: f64,
    pub x: Vec<f64>,
}

impl Coord {
    fn new(range: Range, n: usize) -> Self {
        assert!(n > 0, "n should be positive value.");
        let (x1, x2) = range.tuple();
        let h = (x2 - x1) / (n as f64);
        Coord {
            range,
            n,
            h,
            x: (0..n + 1).map(|i| x1 + h * (i as f64)).collect(),
        }
    }
}

#[derive(Debug)]
pub struct Grid {
    pub depth: u32,
    pub coord: Coord,
    pub val: Vec<f64>,
}

impl Grid {
    fn new(depth: u32, coord: Coord, val: Vec<f64>) -> Self {
        Grid { depth, coord, val }
    }

    fn make_coord(range: Range, depth: u32) -> Coord {
        Coord::new(
            range,
            2_usize.checked_pow(depth).expect(&format!(
                "2 to the {} is beyond the upper bound of usize type.",
                depth
            )),
        )
    }

    pub fn new_zeros(range: Range, depth: u32) -> Self {
        let coord = Grid::make_coord(range, depth);
        let val = vec![0.; coord.n + 1];
        Grid::new(depth, coord, val)
    }

    pub fn new_func(range: Range, depth: u32, func: fn(&f64) -> f64) -> Self {
        let coord = Grid::make_coord(range, depth);
        let val = coord.x.iter().map(func).collect();
        Grid::new(depth, coord, val)
    }

    pub fn rms(&self) -> f64 {
        (self.val.iter().map(|x| x * x).sum::<f64>() / ((self.coord.n + 1) as f64)).sqrt()
    }

    pub fn fine(&self) -> Self {
        let coord = Grid::make_coord(self.coord.range.clone(), self.depth + 1);
        let val = (0..(coord.n + 1))
            .map(|i| {
                if i % 2 == 0 {
                    self.val[i / 2]
                } else {
                    self.val[i / 2] / 2. + self.val[i / 2 + 1] / 2.
                }
            })
            .collect();
        Grid::new(self.depth + 1, coord, val)
    }

    pub fn coarsen(&self) -> Self {
        assert!(
            self.depth > 0,
            "Coarsening is not possible for depth=0 grid."
        );
        let coord = Grid::make_coord(self.coord.range.clone(), self.depth - 1);
        let val = (0..(coord.n + 1))
            .map(|i| {
                if i == 0 {
                    self.val[0]
                } else if i == coord.n {
                    self.val[self.coord.n]
                } else {
                    self.val[2 * i - 1] / 4. + self.val[2 * i] / 2. + self.val[2 * i + 1] / 4.
                }
            })
            .collect();
        Grid::new(self.depth - 1, coord, val)
    }
}
