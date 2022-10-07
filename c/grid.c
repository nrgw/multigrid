#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "grid.h"

Grid *grid_new(double x1, double x2, int n, int num_vars)
{
    Grid *grid = (Grid*)malloc(sizeof(Grid));
    grid->x1 = x1;
    grid->x2 = x2;
    grid->n = n;
    int N = 1 << n;
    grid->N = N;
    double h = (x2 - x1)/(double)N;
    grid->h = h;
    grid->x = (double*)malloc(sizeof(double)*(N + 1));
    for (int i = 0; i <= N; i++) {
        grid->x[i] = x1 + i*h;
    }
    grid->num_vars = num_vars;
    grid->var = (double**)malloc(sizeof(double*)*num_vars);
    for (int i = 0; i < num_vars; i++) {
        grid->var[i] = (double*)malloc(sizeof(double)*(N + 1));
    }
    return grid;
}

void grid_delete(Grid *grid)
{
    for (int i = 0; i < grid->num_vars; i++) {
        free(grid->var[i]);
    }
    free(grid->var);
    free(grid->x);
    free(grid);
}

void grid_set_zero(Grid *grid, int index)
{
    memset(grid->var[index], 0, sizeof(double)*(grid->N + 1));
}

void grid_set_func(Grid *grid, int index, double (*func)(double))
{
    for (int i = 0; i <= grid->N; i++) {
        grid->var[index][i] = func(grid->x[i]);
    }
}

void grid_set_from_coarser(Grid *grid, int index, Grid *coarser, int coarser_index)
{
    // Linear Interpolation
    for (int i = 0; i <= grid->N; i++) {
        if (i % 2 == 0) {
            grid->var[index][i] = coarser->var[coarser_index][i/2];
        } else {
            grid->var[index][i] = coarser->var[coarser_index][i/2]/2. + coarser->var[coarser_index][i/2 + 1]/2.;
        }
    }
}

void grid_set_from_finer(Grid *grid, int index, Grid *finer, int finer_index)
{
    // Full Weighting Restriction
    grid->var[index][0] = finer->var[finer_index][0];
    for (int i = 1; i < grid->N; i++) {
        grid->var[index][i] = finer->var[finer_index][2*i - 1]/4. + finer->var[finer_index][2*i]/2. + finer->var[finer_index][2*i + 1]/4.;
    
    }
    grid->var[index][grid->N] = finer->var[finer_index][finer->N];
}

void grid_add(Grid *grid, int index, Grid *from, int from_index)
{
    for (int i = 0; i <= grid->N; i++) {
        grid->var[index][i] += from->var[from_index][i];
    }
}

void grid_subtract(Grid *grid, int index, Grid *a, int a_index, Grid *b, int b_index)
{
    for (int i = 0; i <= grid->N; i++) {
        grid->var[index][i] = a->var[a_index][i] - b->var[b_index][i];
    }
}

#define SQ(x) ((x)*(x))

double grid_get_rms(Grid *grid, int index)
{
    double sum = 0.;
    for (int i = 0; i <= grid->N; i++) {
        sum += SQ(grid->var[index][i]);
    }
    return sqrt(sum/(double)(1 + grid->N));
}