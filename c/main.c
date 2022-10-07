#include "multigrid.h"
#include <stdio.h>
#include <time.h>
#include <math.h>

double relax_left_func(double *sol, double *src, double s, double h)
{
    return sol[1] - src[0]*h*h/6.;
}

double relax_middle_func(double *sol, double *src, double s, double h, int i)
{
    return (
        + sol[i + 1]*(1. + h/s)
        + sol[i - 1]*(1. - h/s)
        - src[i]*h*h
    )/2.;
}

double relax_right_func(double *sol, double *src, double s, double h)
{
    return 0.;
}

double res_left_func(double *sol, double *src, double s, double h)
{
    return -(relax_left_func(sol, src, s, h) - sol[0])*6./(h*h);
}

double res_middle_func(double *sol, double *src, double s, double h, int i)
{
    return -(relax_middle_func(sol, src, s, h, i) - sol[i])*2./(h*h);
}

double res_right_func(double *sol, double *src, double s, double h)
{
    return 0.;
}

const double rho_c = 1.28e-3;
const double N = 1.;
const double K = 100.;
#define PI 3.14159265358979323846

double sinc(double x) {
    if (x == 0.) {
        return 1.;
    } else {
        return sin(x)/x;
    }
}

double alpha_sq() {
    return (N + 1.)*K*pow(rho_c, 1./N  - 1.)/(4.* PI);
}

double src_func(double s)
{
    if (s < 0.5) {
        double r_s_sq = PI*PI*alpha_sq();
        double rho = rho_c*sinc(PI*s/(1. - s));
        double a = 1. - s;
        return 4.*PI*rho*r_s_sq/(a*a*a*a);
    } else {
        return 0.;
    }
}

double exact_sol_func(double s)
{
    double factor = 0.;
    if (s < 0.5) {
        factor = 1. + sinc(PI*s/(1. - s));
    } else {
        factor = (1. - s)/s;
    }

    return -4.*PI*rho_c*alpha_sq()*factor;
}

int main()
{
    double s1 = 0.;
    double s2 = 1.;
    BVP *bvp = bvp_new(s1, s2, relax_left_func, relax_middle_func, relax_right_func, res_left_func, res_middle_func, res_right_func, src_func, exact_sol_func);

    int n = 16;
    BVPSolver *solver = bvpsolver_new(bvp, n, 4, 1, 4);
    Grid *exact_sol_grid = grid_new_func(s1, s2, n, exact_sol_func);

    const int number_of_iter = 30;
    clock_t start = clock();
    for (int i = 0; i < number_of_iter; i++) {
        bvpsolver_solve(solver);
        printf("%d %e %e\n", i, bvpsolver_residual_rms_normalized(solver), bvpsolver_error_rms(solver, exact_sol_grid));
    }
    printf("averaged time: %f s\n",(double)(clock() - start)/CLOCKS_PER_SEC/(double)number_of_iter);

    grid_delete(exact_sol_grid);
    bvpsolver_delete(solver);
    bvp_delete(bvp);
    return 0;
}