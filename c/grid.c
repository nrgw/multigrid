#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "grid.h"

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