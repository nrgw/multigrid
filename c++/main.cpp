#include "multigrid.hpp"
#include <iostream>
#include <chrono>
#include <memory>

using namespace std;
using namespace chrono;

const double r_s = 8.;
const double rho_c = 1.28e-3;
const double PI = 3.14159265358979323846;

double src_func(double s)
{
    if (s < 0.5)
    {
        double a = 1. - s;
        double b = s / a;
        double rho = rho_c * (1. - b * b);
        return 4. * PI * rho * (r_s * r_s) / (a * a * a * a);
    }
    else
    {
        return 0.;
    }
}

double exact_sol_func(double s)
{
    if (s < 0.5)
    {
        double a = s / (1. - s);
        return -2. * PI * rho_c * (r_s * r_s) * (1. / 2. - a * a / 3. + a * a * a * a / 10.);
    }
    else
    {
        return -8. / 15. * PI * rho_c * (r_s * r_s) * (1. - s) / s;
    }
}

class BVPSolver : public Multigrid
{
private:
public:
    int n;
    double x1, x2;
    Grid *src_grid;
    Grid *sol_grid;
    Grid *res_grid;
    int num_iter1;
    int num_iter2;
    int num_iter3;
    // Aliases
    double h;
    std::vector<double> *x;
    int N;
    std::vector<double> *sol_val;
    std::vector<double> *src_val;
    std::vector<double> *res_val;

    BVPSolver(double x1_, double x2_, int n_, Grid *src_grid_, int num_iter1_, int num_iter2_, int num_iter3_)
    {
        n = n_;
        x1 = x1_;
        x2 = x2_;
        sol_grid = new Grid(x1, x2, n);
        if (src_grid_ == NULL)
        {
            src_grid = new Grid(x1, x2, n, src_func);
        }
        else
        {
            src_grid = src_grid_;
        }
        res_grid = new Grid(x1, x2, n);

        num_iter1 = num_iter1_;
        num_iter2 = num_iter2_;
        num_iter3 = num_iter3_;

        h = sol_grid->h;
        x = &(sol_grid->x);
        N = sol_grid->N;
        sol_val = &(sol_grid->val);
        src_val = &(src_grid->val);
        res_val = &(res_grid->val);
    }
    ~BVPSolver(){};
    void bvpsolver_relax();
    void bvpsolver_get_residual();
    void bvpsolver_multigrid();
    void bvpsolver_solve();
    Grid *bvpsolver_exact_error();
    double bvpsolver_exact_error_rms();

    double relax_left_func(std::vector<double> *sol, std::vector<double> *src, double s, double h) override
    {
        return sol->at(1) - src->at(0) * h * h / 6.;
    }

    double relax_middle_func(std::vector<double> *sol, std::vector<double> *src, double s, double h, int i) override
    {
        return (
                   +sol->at(i + 1) * (1. + h / s) + sol->at(i - 1) * (1. - h / s) - src->at(i) * h * h) /
               2.;
    }

    double relax_right_func(std::vector<double> *sol, std::vector<double> *src, double s, double h) override
    {
        return 0.;
    }

    double res_left_func(std::vector<double> *sol, std::vector<double> *src, double s, double h) override
    {
        return -(relax_left_func(sol, src, s, h) - sol->at(0)) * 6. / (h * h);
    }

    double res_middle_func(std::vector<double> *sol, std::vector<double> *src, double s, double h, int i) override
    {
        return -(relax_middle_func(sol, src, s, h, i) - sol->at(i)) * 2. / (h * h);
    }

    double res_right_func(std::vector<double> *sol, std::vector<double> *src, double s, double h) override
    {
        return 0.;
    }
};

void BVPSolver::bvpsolver_relax()
{
    sol_val->at(0) = relax_left_func(sol_val, src_val, x->at(0), h);
    // Residual at middle points
    for (int i = 1; i < N; i++)
        sol_val->at(i) = relax_middle_func(sol_val, src_val, x->at(i), h, i);
    // Residual at the rightmost point
    sol_val->at(N) = relax_right_func(sol_val, src_val, x->at(N), h);
}
void BVPSolver::bvpsolver_get_residual()
{
    res_val->at(0) = res_left_func(sol_val, src_val, x->at(0), h);
    // Residual at middle points
    for (int i = 1; i < N; i++)
        res_val->at(i) = res_middle_func(sol_val, src_val, x->at(i), h, i);
    // Residual at the rightmost point
    res_val->at(N) = res_right_func(sol_val, src_val, x->at(N), h);
}
void BVPSolver::bvpsolver_multigrid()
{
    for (int j = 0; j < num_iter1; j++)
    {
        bvpsolver_relax();
    }
    if (n != 1)
    {
        for (int j = 0; j < num_iter2; j++)
        {
            bvpsolver_get_residual();
            BVPSolver *coarse_solver = new BVPSolver(x1, x2, n - 1, res_grid->Grid_coarsen(), num_iter1, num_iter2, num_iter3);
            coarse_solver->bvpsolver_solve();
            Grid *error_grid = coarse_solver->sol_grid->Grid_fine();

            for (int i = 0; i <= N; i++)
            {
                sol_val->at(i) += error_grid->val[i];
            }
            delete error_grid;
            delete coarse_solver;
        }

        for (int j = 0; j < num_iter3; j++)
        {
            bvpsolver_relax();
        }
    }
}
void BVPSolver::bvpsolver_solve()
{
    bvpsolver_multigrid();
}
Grid *BVPSolver::bvpsolver_exact_error()
{
    if (exact_sol_func == NULL)
    {
        std::cout << "exact_sol_func is not specified.\n";
        exit(1);
    }

    Grid *error_grid = new Grid(x1, x2, n, exact_sol_func);
    for (int i = 0; i <= N; i++)
    {
        error_grid->val[i] = error_grid->val[i] - sol_val->at(i);
    }

    return error_grid;
}
double BVPSolver::bvpsolver_exact_error_rms()
{
    Grid *error_grid = bvpsolver_exact_error();
    double error_rms = error_grid->grid_rms();
    delete error_grid;
    return error_rms;
}

template <typename... Args>
std::string string_format(const std::string &format, Args... args)
{
    size_t size = snprintf(nullptr, 0, format.c_str(), args...) + 1; // Extra space for '\0'
    if (size <= 0)
    {
        throw std::runtime_error("Error during formatting.");
    }
    std::unique_ptr<char[]> buf(new char[size]);
    snprintf(buf.get(), size, format.c_str(), args...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

int main()
{
    double s1 = 0.;
    double s2 = 1.;

    const int number_of_iter = 30;

    BVPSolver solver = BVPSolver(s1, s2, 16, NULL, 4, 1, 4);

    system_clock::time_point start = system_clock::now();
    for (int i = 0; i < number_of_iter; i++)
    {
        solver.bvpsolver_solve();
        cout << string_format("%d %e %e", i, solver.res_grid->grid_rms(), solver.bvpsolver_exact_error_rms()) << endl;
    }
    system_clock::time_point end = system_clock::now();
    seconds m_time = duration_cast<seconds>(end - start);

    cout << string_format("averaged time: %f s", m_time.count() / (double)number_of_iter) << endl;
    return 0;
}
