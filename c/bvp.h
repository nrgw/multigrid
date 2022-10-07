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

typedef enum _GridKind {SRC, SOL, RES, ERR, NUM_GK} GridKind;

typedef struct {
    BVP *bvp;
    int n;
    Grid **grid[NUM_GK];
    int num_iter1;
    int num_iter2;
    int num_iter3;
} BVPSolver;

BVPSolver *bvpsolver_new(BVP *bvp, int num_level, int num_iter1, int num_iter2, int num_iter3);
void bvpsolver_delete(BVPSolver *solver);
void bvpsolver_solve(BVPSolver *solver);
double bvpsolver_residual_rms(BVPSolver *solver);
double bvpsolver_residual_rms_normalized(BVPSolver *solver);
double bvpsolver_error_rms(BVPSolver *solver, Grid *solution);

#endif