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
    memset(grid->val, 0, sizeof(double)*(grid->N + 1));
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

double grid_rms(Grid *grid)
{
    double sum = 0.;
    for (int i = 0; i <= grid->N; i++) {
        sum += grid->val[i]*grid->val[i];
    }
    return sqrt(sum/(double)(1 + grid->N));
}

Grid *grid_fine(Grid *grid)
{
    Grid *finer = grid_new(grid->x1, grid->x2, grid->n + 1);
    // Linear Interpolation
    for (int i = 0; i <= finer->N; i++) {
        if (i % 2 == 0) {
            finer->val[i] = grid->val[i/2];
        } else {
            finer->val[i] = grid->val[i/2]/2. + grid->val[i/2 + 1]/2.;
        }
    }
    return finer;
}

Grid *grid_coarsen(Grid *grid)
{
    Grid *coarser = grid_new(grid->x1, grid->x2, grid->n - 1);
    // Full Weighting Restriction
    coarser->val[0] = grid->val[0];
    for (int i = 1; i < coarser->N; i++) {
        coarser->val[i] = grid->val[2*i - 1]/4. + grid->val[2*i]/2. + grid->val[2*i + 1]/4.;
    
    }
    coarser->val[coarser->N] = grid->val[grid->N];
    return coarser;
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

BVPSolver *bvpsolver_new(BVP *bvp, int n, Grid *src_grid, int num_iter1, int num_iter2, int num_iter3)
{
    BVPSolver *solver = (BVPSolver*)malloc(sizeof(BVPSolver));
    solver->bvp = bvp;
    solver->n = n;
    solver->sol_grid = grid_new_zeros(bvp->x1, bvp->x2, n);
    if (src_grid == NULL) {
        solver->src_grid = grid_new_func(bvp->x1, bvp->x2, n, bvp->src_func);
    } else {
        solver->src_grid = src_grid;
    }
    solver->res_grid = grid_new(bvp->x1, bvp->x2, n);
    
    solver->num_iter1 = num_iter1;
    solver->num_iter2 = num_iter2;
    solver->num_iter3 = num_iter3;

    // Aliases
    solver->h = solver->sol_grid->h;
    solver->x = solver->sol_grid->x;
    solver->N = solver->sol_grid->N;
    solver->sol_val = solver->sol_grid->val;
    solver->src_val = solver->src_grid->val;
    solver->res_val = solver->res_grid->val;

    return solver;
}

void bvpsolver_delete(BVPSolver *solver)
{
    grid_delete(solver->src_grid);
    grid_delete(solver->sol_grid);
    grid_delete(solver->res_grid);
    free(solver);
}

void bvpsolver_relax(BVPSolver *solver)
{
    // Aliases
    double h = solver->h;
    double *x = solver->x;
    int N = solver->N;
    double *sol_val = solver->sol_val;
    double *src_val = solver->src_val;

    // Red Sweep
    sol_val[0] = solver->bvp->relax_left_func(sol_val, src_val, x[0], h);
    for (int i = 2; i < N; i += 2)
        sol_val[i] = solver->bvp->relax_middle_func(sol_val, src_val, x[i], h, i);
    sol_val[N] = solver->bvp->relax_right_func(sol_val, src_val, x[N], h);

    // Black Sweep
    for (int i = 1; i < N; i += 2)
        sol_val[i] = solver->bvp->relax_middle_func(sol_val, src_val, x[i], h, i);
}

void bvpsolver_get_residual(BVPSolver *solver)
{
    // Aliases
    double h = solver->h;
    double *x = solver->x;
    int N = solver->N;
    double *sol_val = solver->sol_val;
    double *src_val = solver->src_val;
    double *res_val = solver->res_val;


    // Residual at the leftmost point
    res_val[0] = solver->bvp->res_left_func(sol_val, src_val, x[0], h);
    // Residual at middle points
    for (int i = 1; i < N; i++)
        res_val[i] = solver->bvp->res_middle_func(sol_val, src_val, x[i], h, i);
    // Residual at the rightmost point
    res_val[N] = solver->bvp->res_left_func(sol_val, src_val, x[N], h);
}

void bvpsolver_multigrid(BVPSolver *solver)
{
    // Pre-Smoothing
    for (int j = 0; j < solver->num_iter1; j++)
        bvpsolver_relax(solver);
    if (solver->n != 1) {
        // Error Correction using Coarse Grid
        for (int j = 0; j < solver->num_iter2; j++) {
            bvpsolver_get_residual(solver);
            BVPSolver *coarse_solver = bvpsolver_new(solver->bvp, solver->n - 1, grid_coarsen(solver->res_grid), solver->num_iter1, solver->num_iter2, solver->num_iter3);
            bvpsolver_solve(coarse_solver);
            Grid *error_grid = grid_fine(coarse_solver->sol_grid);
            for (int i = 0; i <= solver->N; i++) {
                solver->sol_val[i] += error_grid->val[i];
            }
            grid_delete(error_grid);
            bvpsolver_delete(coarse_solver);
        }

        // Post-Smoothing
        for (int j = 0; j < solver->num_iter3; j++)
            bvpsolver_relax(solver);
    }
}

void bvpsolver_solve(BVPSolver *solver)
{
    bvpsolver_multigrid(solver);
}

Grid *bvpsolver_exact_error(BVPSolver *solver)
{
    if (solver->bvp->exact_sol_func == NULL) {
        fprintf(stderr, "exact_sol_func is not specified.\n");
        exit(1);
    }

    Grid *error_grid = grid_new_func(solver->bvp->x1, solver->bvp->x2, solver->n, solver->bvp->exact_sol_func);
    for (int i = 0; i <= solver->N; i++) {
        error_grid->val[i] -= solver->sol_val[i];
    }

    return error_grid;
}

double bvpsolver_exact_error_rms(BVPSolver *solver)
{
    Grid *error_grid = bvpsolver_exact_error(solver);
    double error_rms = grid_rms(error_grid);
    grid_delete(error_grid);
    return error_rms;
}