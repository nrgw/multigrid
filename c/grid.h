#ifndef GRID_H
#define GRID_H

typedef struct {
    double x1;
    double x2;
    int n;
    int N;
    double h;
    double *x;
    int num_vars;
    double **var;
} Grid;

Grid *grid_new(double x1, double x2, int n, int num_var);
void grid_delete(Grid *grid);
void grid_set_zero(Grid *grid, int index);
void grid_set_func(Grid *grid, int index, double (*func)(double));
void grid_set_from_coarser(Grid *grid, int index, Grid *coarser, int coarser_index);
void grid_set_from_finer(Grid *grid, int index, Grid *finer, int finer_index);
void grid_add(Grid *grid, int index, Grid *from, int from_index);
void grid_subtract(Grid *grid, int index, Grid *a, int a_index, Grid *b, int b_index);
double grid_get_rms(Grid *grid, int index);

#endif