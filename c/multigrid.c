#include "multigrid.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

Grid *grid_new(double x1, double x2, int n)
{
    Grid *grid = (Grid*)malloc(sizeof(Grid));
    grid->x1 = x1;
    grid->x2 = x2;
    grid->n = n;
    int N = 1<<n;
    grid->N = N;
    double h = (x2 - x1)/(double)N;
    grid->h = h;
    grid->x = (double*)malloc(sizeof(double)*(N + 1));
    for (int i = 0; i <= N; i++) {
        grid->x[i] = x1 + i*h;
    }
    grid->val = malloc(sizeof(double)*(N + 1));
    return grid;
}

Grid *grid_new_zeros(double x1, double x2, int n)
{
    Grid *grid = grid_new(x1, x2, n);
    grid_set_zero(grid);
    return grid;
}

Grid *grid_new_func(double x1, double x2, int n, double (*func)(double))
{
    Grid *grid = grid_new(x1, x2, n);
    for (int i = 0; i <= grid->N; i++) {
        grid->val[i] = func(grid->x[i]);
    }
    return grid;
}

void grid_delete(Grid *grid)
{
    free(grid->x);
    free(grid->val);
    free(grid);
}

void grid_set_zero(Grid *grid)
{
    memset(grid->val, 0, sizeof(double)*(grid->N + 1));
}

double grid_rms(Grid *grid)
{
    double sum = 0.;
    for (int i = 0; i <= grid->N; i++) {
        sum += grid->val[i]*grid->val[i];
    }
    return sqrt(sum/(double)(1 + grid->N));
}

void grid_fine(Grid *grid, Grid *target)
{
    // Linear Interpolation
    for (int i = 0; i <= target->N; i++) {
        if (i % 2 == 0) {
            target->val[i] = grid->val[i/2];
        } else {
            target->val[i] = grid->val[i/2]/2. + grid->val[i/2 + 1]/2.;
        }
    }
}

void grid_coarsen(Grid *grid, Grid *target)
{
    // Full Weighting Restriction
    target->val[0] = grid->val[0];
    for (int i = 1; i < target->N; i++) {
        target->val[i] = grid->val[2*i - 1]/4. + grid->val[2*i]/2. + grid->val[2*i + 1]/4.;
    
    }
    target->val[target->N] = grid->val[grid->N];
}

void grid_subtract(Grid *a, Grid *b, Grid *target)
{
    for (int i = 0; i <= target->N; i++) {
        target->val[i] = a->val[i] - b->val[i];
    }
}

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

BVPSolver *bvpsolver_new(BVP *bvp, int num_level, int num_iter1, int num_iter2, int num_iter3)
{
    BVPSolver *solver = (BVPSolver*)malloc(sizeof(BVPSolver));
    solver->bvp = bvp;
    solver->n = num_level;

    for (GridKind i = 0; i < NUM_GK; i++) {
        solver->grid[i] = (Grid**)malloc(sizeof(Grid*)*num_level);
        for (int j = 0; j < num_level - 1; j++) {
            solver->grid[i][j] = grid_new(bvp->x1, bvp->x2, j + 1);
        }
        switch (i) {
            case SRC:
                solver->grid[i][num_level - 1] = grid_new_func(bvp->x1, bvp->x2, num_level, bvp->src_func);
                break;
            case SOL:
                solver->grid[i][num_level - 1] = grid_new_zeros(bvp->x1, bvp->x2, num_level);
                break;
            case RES:
            case ERR:
                solver->grid[i][num_level - 1] = grid_new(bvp->x1, bvp->x2, num_level);
                break;
        }
    }

    solver->num_iter1 = num_iter1;
    solver->num_iter2 = num_iter2;
    solver->num_iter3 = num_iter3;

    return solver;
}

void bvpsolver_delete(BVPSolver *solver)
{
    for (GridKind i = 0; i < NUM_GK; i++) {
        for (int j = 0; j < solver->n; j++) {
            grid_delete(solver->grid[i][j]);
        }
        free(solver->grid[i]);
    }
    free(solver);
}

void bvpsolver_relax(BVPSolver *solver, int level)
{
    // Aliases
    double h = solver->grid[0][level]->h;
    double *x = solver->grid[0][level]->x;
    int N = solver->grid[0][level]->N;
    double *src_val = solver->grid[SRC][level]->val;
    double *sol_val = solver->grid[SOL][level]->val;

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
    double h = solver->grid[0][level]->h;
    double *x = solver->grid[0][level]->x;
    int N = solver->grid[0][level]->N;
    double *src_val = solver->grid[SRC][level]->val;
    double *sol_val = solver->grid[SOL][level]->val;
    double *res_val = solver->grid[RES][level]->val;

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
            grid_coarsen(solver->grid[RES][level], solver->grid[SRC][level - 1]);
            grid_set_zero(solver->grid[SOL][level - 1]);
            bvpsolver_multigrid(solver, level - 1);
            grid_fine(solver->grid[SOL][level - 1], solver->grid[ERR][level]);
            for (int i = 0; i <= solver->grid[SOL][level]->N; i++) {
                solver->grid[SOL][level]->val[i] += solver->grid[ERR][level]->val[i];
            }
        }

        // Post-Smoothing
        for (int j = 0; j < solver->num_iter3; j++)
            bvpsolver_relax(solver, level);
    }
}

void bvpsolver_get_error(BVPSolver *solver, Grid *solution)
{
    grid_subtract(solution, solver->grid[SOL][solver->n - 1], solver->grid[ERR][solver->n - 1]);
}

// public functions

void bvpsolver_solve(BVPSolver *solver)
{
    bvpsolver_multigrid(solver, solver->n - 1);
}

double bvpsolver_residual_rms(BVPSolver *solver)
{
    bvpsolver_get_residual(solver, solver->n - 1);
    return grid_rms(solver->grid[RES][solver->n - 1]);
}

double bvpsolver_residual_rms_normalized(BVPSolver *solver)
{
    double h = solver->grid[SOL][solver->n - 1]->h;
    return bvpsolver_residual_rms(solver)*h*h;
}

double bvpsolver_error_rms(BVPSolver *solver, Grid *solution)
{
    bvpsolver_get_error(solver, solution);
    return grid_rms(solver->grid[ERR][solver->n - 1]);
}