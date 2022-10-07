#include <stdlib.h>
#include "bvp.h"

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
)
{
    BVP *bvp = (BVP*)malloc(sizeof(BVP));
    bvp->x1 = x1;
    bvp->x2 = x2;
    bvp->relax_left_func = relax_left_func;
    bvp->relax_middle_func = relax_middle_func;
    bvp->relax_right_func = relax_right_func;
    bvp->res_left_func = res_left_func;
    bvp->res_middle_func = res_middle_func;
    bvp->res_right_func = res_right_func;
    bvp->src_func = src_func;
    bvp->exact_sol_func = exact_sol_func;
    return bvp;
}

void bvp_delete(BVP *bvp)
{
    free(bvp);
}

typedef enum _GridKind {SOL, SRC, RES, ERR, NUM_GK} GridKind;

BVPSolver *bvpsolver_new(BVP *bvp, int num_level, int num_iter1, int num_iter2, int num_iter3)
{
    BVPSolver *solver = (BVPSolver*)malloc(sizeof(BVPSolver));
    solver->bvp = bvp;
    solver->num_level = num_level;
    solver->subgrid = (Grid**)malloc(sizeof(Grid*)*num_level);
    for (int i = 0; i < num_level; i++) {
        solver->subgrid[i] = grid_new(bvp->x1, bvp->x2, i + 1, NUM_GK);
    }
    solver->grid = solver->subgrid[num_level - 1];
    grid_set_zero(solver->grid, SOL);
    grid_set_func(solver->grid, SRC, bvp->src_func);

    solver->num_iter1 = num_iter1;
    solver->num_iter2 = num_iter2;
    solver->num_iter3 = num_iter3;

    return solver;
}

void bvpsolver_delete(BVPSolver *solver)
{
    for (int i = 0; i < solver->num_level; i++) {
        grid_delete(solver->subgrid[i]);
    }
    free(solver->subgrid);
    free(solver);
}

void bvpsolver_relax(BVPSolver *solver, int level)
{
    // Aliases
    double h = solver->subgrid[level]->h;
    double *x = solver->subgrid[level]->x;
    int N = solver->subgrid[level]->N;
    double *sol_val = solver->subgrid[level]->var[SOL];
    double *src_val = solver->subgrid[level]->var[SRC];

    // Red Sweep
    sol_val[0] = solver->bvp->relax_left_func(sol_val, src_val, x[0], h);
    for (int i = 2; i < N; i += 2)
        sol_val[i] = solver->bvp->relax_middle_func(sol_val, src_val, x[i], h, i);
    sol_val[N] = solver->bvp->relax_right_func(sol_val, src_val, x[N], h);

    // Black Sweep
    for (int i = 1; i < N; i += 2)
        sol_val[i] = solver->bvp->relax_middle_func(sol_val, src_val, x[i], h, i);
}

void bvpsolver_get_residual(BVPSolver *solver, int level)
{
    // Aliases
    double h = solver->subgrid[level]->h;
    double *x = solver->subgrid[level]->x;
    int N = solver->subgrid[level]->N;
    double *sol_val = solver->subgrid[level]->var[SOL];
    double *src_val = solver->subgrid[level]->var[SRC];
    double *res_val = solver->subgrid[level]->var[RES];

    // Residual at the leftmost point
    res_val[0] = solver->bvp->res_left_func(sol_val, src_val, x[0], h);
    // Residual at middle points
    for (int i = 1; i < N; i++)
        res_val[i] = solver->bvp->res_middle_func(sol_val, src_val, x[i], h, i);
    // Residual at the rightmost point
    res_val[N] = solver->bvp->res_left_func(sol_val, src_val, x[N], h);
}

void bvpsolver_multigrid(BVPSolver *solver, int level)
{
    // Pre-Smoothing
    for (int j = 0; j < solver->num_iter1; j++)
        bvpsolver_relax(solver, level);
    if (level > 0) {
        // Error Correction using Coarse Grid
        for (int j = 0; j < solver->num_iter2; j++) {
            bvpsolver_get_residual(solver, level);
            grid_coarsen(solver->subgrid[level], RES, solver->subgrid[level - 1], SRC);
            grid_set_zero(solver->subgrid[level - 1], SOL);
            bvpsolver_multigrid(solver, level - 1);
            grid_fine(solver->subgrid[level - 1], SOL, solver->subgrid[level], ERR);
            grid_add(solver->subgrid[level], SOL, solver->subgrid[level], ERR);
        }

        // Post-Smoothing
        for (int j = 0; j < solver->num_iter3; j++)
            bvpsolver_relax(solver, level);
    }
}

void bvpsolver_get_error(BVPSolver *solver, Grid *solution)
{
    grid_subtract(solver->grid, ERR, solution, 0, solver->grid, SOL);
}

// public functions

void bvpsolver_solve(BVPSolver *solver)
{
    bvpsolver_multigrid(solver, solver->num_level - 1);
}

double bvpsolver_residual_rms(BVPSolver *solver)
{
    bvpsolver_get_residual(solver, solver->num_level - 1);
    return grid_rms(solver->grid, RES);
}

#define SQ(x) ((x)*(x))

double bvpsolver_residual_rms_normalized(BVPSolver *solver)
{
    return bvpsolver_residual_rms(solver)*SQ(solver->grid->h);
}

double bvpsolver_error_rms(BVPSolver *solver, Grid *solution)
{
    bvpsolver_get_error(solver, solution);
    return grid_rms(solver->grid, ERR);
}