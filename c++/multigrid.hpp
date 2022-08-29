#ifndef __MULTIGRID_H__
#define __MULTIGRID_H__
#include <vector>
#include <iostream>

class Grid
{
public:
    double x1;
    double x2;
    int n;
    int N;
    double h;
    std::vector<double> x;
    std::vector<double> val;

    Grid(double x1_, double x2_, int n_) : x1{x1_}, x2{x2_}, n{n_}
    {
        init();
    }
    Grid(double x1_, double x2_, int n_, double (*func)(double)) : x1{x1_}, x2{x2_}, n{n_}
    {
        init();
        for (auto i = 0; i <= N; i++)
        {
            val[i] = func(x[i]);
        }
    }
    void init()
    {
        N = 1 << n;
        h = (x2 - x1) / (double)N;
        for (auto i = 0; i <= N; i++)
        {
            x.push_back(x1 + i * h);
        }
        val.assign(N + 1, 0);
    }

    Grid *Grid_fine();
    Grid *Grid_coarsen();

    double grid_rms();
    void show();
};

class Multigrid
{
public:
    virtual double relax_left_func(std::vector<double> *, std::vector<double> *, double, double) = 0;
    virtual double relax_middle_func(std::vector<double> *, std::vector<double> *, double, double, int) = 0;
    virtual double relax_right_func(std::vector<double> *, std::vector<double> *, double, double) = 0;
    virtual double res_left_func(std::vector<double> *, std::vector<double> *, double, double) = 0;
    virtual double res_middle_func(std::vector<double> *, std::vector<double> *, double, double, int) = 0;
    virtual double res_right_func(std::vector<double> *, std::vector<double> *, double, double) = 0;
    virtual ~Multigrid() = default;
};

#endif