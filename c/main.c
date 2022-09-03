#include "multigrid.h"
#include <stdio.h>
#include <time.h>

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

const double r_s = 8.;
const double rho_c = 1.28e-3;
#define PI 3.14159265358979323846

double src_func(double s)
{
    if (s < 0.5) {
        double a = 1. - s;
        double b = s / a;
        double rho = rho_c*(1. - b*b);
        return 4.*PI*rho*(r_s*r_s)/(a*a*a*a);
    } else {
        return 0.;
    }
}

double exact_sol_func(double s)
{
    if (s < 0.5) {
        double a = s / (1. - s);
        return -2.*PI*rho_c*(r_s*r_s)*(1./2. - a*a/3. + a*a*a*a/10.);
    } else {
        return -8./15.*PI*rho_c*(r_s*r_s)*(1. - s)/s;
    }
}

int main()
{
    double s1 = 0.;
    double s2 = 1.;
    BVP *bvp = bvp_new(s1, s2, relax_left_func, relax_middle_func, relax_right_func, res_left_func, res_middle_func, res_right_func, src_func, exact_sol_func);
    BVPSolver *solver = bvpsolver_new(bvp, 16, NULL, 4, 1, 4);

    const int number_of_iter = 30;
    clock_t start = clock();
    for (int i = 0; i < number_of_iter; i++) {
        bvpsolver_solve(solver);
        printf("%d %e %e\n", i, grid_rms(solver->res_grid)*solver->res_grid->h, bvpsolver_exact_error_rms(solver));
    }
    printf("averaged time: %f s\n",(double)(clock() - start)/CLOCKS_PER_SEC/(double)number_of_iter);

    bvpsolver_delete(solver);
    bvp_delete(bvp);
    return 0;
}
