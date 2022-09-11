import numpy as np

def relax_left(sol, src, s, h):
    return sol[1] - src[0]*h**2/6

def relax_middle(sol, src, s, h, i):
    return (
        + sol[i + 1]*(1 + h/s)
        + sol[i - 1]*(1 - h/s)
        - src[i]*h**2
    )/2

def relax_right(sol, src, s, h):
    return 0

def residual_left(sol, src, s, h):
    return -(relax_left(sol, src, s, h) - sol[0])*6/h**2

def residual_middle(sol, src, s, h, i):
    return -(relax_middle(sol, src, s, h, i) - sol[i])*2/h**2

def residual_right(sol, src, s, h):
    return 0

r_s = 8
rho_c = 1.28e-3

def src(s):
    if s < 0.5:
        rho = rho_c*(1 - (s / (1 - s))**2)
        return 4*np.pi*rho*r_s**2*(1 - s)**(-4)
    else:
        return 0

def exact_sol(s):
    if s < 0.5:
        a = s/(1 - s)
        return -2*np.pi*rho_c*r_s**2*(1/2 - a**2/3 + a**4/10)
    else:
        return -8/15*np.pi*rho_c*r_s**2*(1 - s)/s

from multigrid import BVP, BVPSolver, BVPSolver_vec

import time

def run_bpv_solver(BVPSolver=BVPSolver):
    s1 = 0
    s2 = 1
    bvp = BVP((s1, s2),
              relax_left, relax_middle, relax_right,
              residual_left, residual_middle, residual_right,
              src, exact_sol_func=exact_sol)

    n = 16
    solver = BVPSolver(bvp, n, num_iter=(4,1,4))

    number_of_iter = 10
    start = time.time()
    for i in range(number_of_iter):
        solver.solve()
        res_rms = solver.residual().rms()
        print(i, res_rms*solver.sol_grid.h**2, solver.exact_error().rms())
    print("time : {} s".format((time.time() - start)/number_of_iter))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--vec', dest='solver', action='store_const',
                        const=BVPSolver_vec, default=BVPSolver,
                        help='use vectorized version')

    args = parser.parse_args()

    run_bpv_solver(args.solver)
