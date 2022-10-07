import time
import numpy as np
from multigrid import BVP, Grid


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

rho_c = 1.28e-3
N = 1
K = 100
alpha_sq = (N + 1)*K*rho_c**(1/N  - 1)/(4*np.pi)
r_s_sq = np.pi**2*alpha_sq

def src(s):
    if s < 0.5:
        rho = rho_c*np.sinc(s/(1 - s))
        return 4*np.pi*rho*r_s_sq/(1 - s)**4
    else:
        return 0

def exact_sol(s):
    if s < 0.5:
        factor = 1 + np.sinc(s/(1 - s))
    else:
        factor = (1 - s)/s
    return -4*np.pi*rho_c*alpha_sq*factor

s1 = 0
s2 = 1
bvp = BVP((s1, s2),
          relax_left, relax_middle, relax_right,
          residual_left, residual_middle, residual_right,
          src, exact_sol_func=exact_sol)


def run_bpv_solver(bvp, BVPSolver, depth=16, niter=10):

    n = depth
    Grid.preallocate_memory(n, bvp.domain)
    solver = BVPSolver(bvp, n, num_iter=(4,1,4))

    number_of_iter = niter
    times = []
    for i in range(number_of_iter):
        start = time.time()
        solver.solve()
        res_rms = solver.residual().rms()
        times.append(time.time() - start)
        print(i, res_rms*solver.sol_grid.h**2, solver.exact_error().rms())

    print("time : {} s".format(np.median(times)))


if __name__ == '__main__':
    import argparse

    from multigrid import BVPSolver
    from multigrid_numpy import BVPSolver_numpy
    try:
        import jax
    except ImportError:
        print("Importing jax failed. The jax-solver is disabled.")
        jax = None
    try:
        import numba
    except ImportError:
        print("Importing numba failed. The jax-solver is disabled.")
        numba = None

    solvers = dict(default=BVPSolver,
                   numpy=BVPSolver_numpy)

    if jax is not None:
        from multigrid_jax import BVPSolver_jax
        solvers["jax"] = BVPSolver_jax

    if numba is not None:
        from multigrid_numba import BVPSolver_numba
        solvers["numba"] = BVPSolver_numba

    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', nargs=1, help='solver to use',
                        choices=solvers.keys(), type=str,
                        dest="solver", default=["default"])

    parser.add_argument('--depth', help='depth of grid',
                        type=int,
                        dest="depth", default=16)

    parser.add_argument('--niter', help='nubner of iter',
                        type=int,
                        dest="niter", default=10)

    args = parser.parse_args()

    solver_name = args.solver[0]
    solver = solvers[solver_name]

    print(f"solver:{solver_name}, depth={args.depth}, niter={args.niter}")
    run_bpv_solver(bvp, solver, depth=args.depth, niter=args.niter)
