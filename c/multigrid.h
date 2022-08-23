typedef struct {
    double x1;
    double x2;
    int n;
    int N;
    double h;
    double *x;
    double *val;
} Grid;

Grid *grid_new(double x1, double x2, int n);
Grid *grid_new_zeros(double x1, double x2, int n);
Grid *grid_new_func(double x1, double x2, int n, double (*func)(double));
void grid_delete(Grid *grid);
double grid_rms(Grid *grid);
Grid *grid_fine(Grid *grid);
Grid *grid_coarsen(Grid *grid);

typedef struct {
    double x1;
    double x2;
    double (*relax_left_func)(double*, double*, double, double);
    double (*relax_middle_func)(double*, double*, double, double, int);
    double (*relax_right_func)(double*, double*, double, double);
    double (*res_left_func)(double*, double*, double, double);
    double (*res_middle_func)(double*, double*, double, double, int);
    double (*res_right_func)(double*, double*, double, double);
    double (*src_func)(double);
    double (*exact_sol_func)(double);  
} BVP;

BVP *bvp_new(
    double x1, double x2,
    double (*relax_left_func)(double*, double*, double, double),
    double (*relax_middle_func)(double*, double*, double, double, int),
    double (*relax_right_func)(double*, double*, double, double),
    double (*res_left_func)(double*, double*, double, double),
    double (*res_middle_func)(double*, double*, double, double, int),
    double (*res_right_func)(double*, double*, double, double),
    double (*src_func)(double),
    double (*exact_sol_func)(double)
);

void bvp_delete(BVP *bvp);

typedef struct {
    BVP *bvp;
    int n;
    Grid *src_grid;
    Grid *sol_grid;
    Grid *res_grid;
    int num_iter1;
    int num_iter2;
    int num_iter3;
    // Aliases
    double h;
    double *x;
    int N;
    double *sol_val;
    double *src_val;
    double *res_val;
} BVPSolver;

BVPSolver *bvpsolver_new(BVP *bvp, int n, Grid *src_grid, int num_iter1, int num_iter2, int num_iter3);
void bvpsolver_delete(BVPSolver *solver);
void bvpsolver_relax(BVPSolver *solver);
void bvpsolver_get_residual(BVPSolver *solver);
void bvpsolver_multigrid(BVPSolver *solver);
void bvpsolver_solve(BVPSolver *solver);
Grid *bvpsolver_exact_error(BVPSolver *solver);
double bvpsolver_exact_error_rms(BVPSolver *solver);