#ifndef BVP_H
#define BVP_H

#include "grid.h"

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
    int num_levels;
    Grid *top_grid;
    Grid **grid;
    int num_iters[3];
} BVPSolver;

BVPSolver *bvpsolver_new(BVP *bvp, int num_levels, int num_iters_1, int num_iters_2, int num_iters_3);
void bvpsolver_delete(BVPSolver *solver);
void bvpsolver_solve(BVPSolver *solver);
double bvpsolver_residual_rms(BVPSolver *solver);
double bvpsolver_residual_rms_normalized(BVPSolver *solver);
double bvpsolver_error_rms(BVPSolver *solver, Grid *solution);

#endif