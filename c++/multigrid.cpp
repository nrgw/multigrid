#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <iostream>
#include "multigrid.hpp"

double Grid::grid_rms()
{
    double sum = 0.;
    for (int i = 0; i <= N; i++)
    {
        sum += val[i] * val[i];
    }

    return std::sqrt(sum / (double)(N + 1));
}

Grid *Grid::Grid_fine()
{
    Grid *finer = new Grid(x1, x2, n + 1);

    for (int i = 0; i <= finer->N; i++)
    {
        if (i % 2 == 0)
        {
            finer->val[i] = val[i / 2];
        }
        else
        {
            finer->val[i] = val[i / 2] / 2. + val[i / 2 + 1] / 2.;
        }
    }
    return finer;
}

Grid *Grid::Grid_coarsen()
{
    Grid *coarser = new Grid(x1, x2, n - 1);

    coarser->val[0] = val[0];
    for (int i = 0; i < coarser->N; i++)
    {
        coarser->val[i] = val[2 * i - 1] / 4. + val[2 * i] / 2. + val[2 * i + 1] / 4.;
    }
    coarser->val[coarser->N] = val[N];
    return coarser;
}

void Grid::show()
{
    for (double v : x)
    {
        std::cout << v << " ";
    }
}