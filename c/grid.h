#ifndef GRID_H
#define GRID_H

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
void grid_set_zero(Grid *grid);
double grid_rms(Grid *grid);
void grid_fine(Grid *grid, Grid *target);
void grid_coarsen(Grid *grid, Grid *target);
void grid_subtract(Grid *a, Grid *b, Grid *target);

#endif