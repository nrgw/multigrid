from scipy.interpolate import interp1d
import numpy as np

class Grid:
    def __init__(self, domain, n, val=None, func=None, shape=None):
        self.domain = domain
        self.n = n

        N = 2**n
        x1 = domain[0]
        x2 = domain[1]

        self.N = N
        self.x = np.linspace(x1, x2, num=N+1)
        self.h = (x2 - x1)/N
        
        if val is not None:
            if val.shape[0] != N + 1:
                raise Exception('Shape of val is invalid.')
            self.val = val
        elif func is not None:
            self.val = np.array([func(xi) for xi in self.x])
        elif shape is not None:
            self.val = np.zeros((self.N + 1,) + shape)
        else:
            raise Exception('One of val, func, and shape should be specified.')
            
    def rms(self):
        return np.linalg.norm(self.val, ord=2, axis=0)/np.sqrt(self.val.shape[0])

    def fine(self, kind='linear'):
        finer = Grid(self.domain, self.n + 1, shape=self.val.shape[1:])
        finer.val = interp1d(self.x, self.val, kind=kind, axis=0)(finer.x)
        return finer

    def coarsen(self):
        val = self.val
        coarser = Grid(self.domain, self.n - 1, shape=val.shape[1:])
        # Full Weighting Restriction
        coarser.val[0] = val[0]
        coarser.val[-1] = val[-1]
        for i in range(1, coarser.N):
            coarser.val[i] = 1/4*val[2*i - 1] + 1/2*val[2*i] + 1/4*val[2*i + 1]
        return coarser

class BVP:
    def __init__(self, domain, relax_left_func, relax_middle_func, relax_right_func, res_left_func, res_middle_func, res_right_func, src_func, exact_sol_func=None):
        if not domain[0] < domain[1]:
            raise Exception('Invalid domain.')
        self.domain = domain
        self.relax_left_func = relax_left_func
        self.relax_middle_func = relax_middle_func
        self.relax_right_func = relax_right_func
        self.res_left_func = res_left_func
        self.res_middle_func = res_middle_func
        self.res_right_func = res_right_func
        self.src_func = src_func
        self.exact_sol_func = exact_sol_func

class BVPSolver:
    def __init__(self, bvp, n, src_grid=None, num_iter=(1,1,1)):
        if not isinstance(n, int) and n > 0:
            raise Exception('Invalid n.')

        self.bvp = bvp
        self.n = n

        self.src_grid = Grid(bvp.domain, n, func=bvp.src_func) if src_grid is None else src_grid
        self.sol_grid = Grid(bvp.domain, n, shape=self.src_grid.val.shape[1:])
        
        self.num_iter = num_iter

        self.h = self.sol_grid.h
        self.x = self.sol_grid.x
        self.N = self.sol_grid.N
        self.sol_val = self.sol_grid.val
        self.src_val = self.src_grid.val

    def relax(self):
        h = self.h
        x = self.x
        N = self.N
        sol_val = self.sol_val
        src_val = self.src_val

        # Relaxation at the leftmost point
        sol_val[0] = self.bvp.relax_left_func(sol_val, src_val, x[0], h)
        # Relaxation at middle points
        for i in range(1, N):
            sol_val[i] = self.bvp.relax_middle_func(sol_val, src_val, x[i], h, i)
        # Relaxation at the rightmost point
        sol_val[N] = self.bvp.relax_right_func(sol_val, src_val, x[N], h)

    def residual(self):
        h = self.h
        x = self.x
        N = self.N
        sol_val = self.sol_val
        src_val = self.src_val

        res = np.zeros(sol_val.shape)
        # Residual at the leftmost point
        res[0] = self.bvp.res_left_func(sol_val, src_val, x[0], h)
        # Residual at middle points
        for i in range(1, N):
            res[i] = self.bvp.res_middle_func(sol_val, src_val, x[i], h, i)
        # Residual at the rightmost point
        res[N] = self.bvp.res_left_func(sol_val, src_val, x[N], h)

        return Grid(self.bvp.domain, self.n, val=res)

    def multigrid(self):
        n = self.n

        # Pre-Smoothing
        for i in range(self.num_iter[0]):
            self.relax()
        if n != 1:
            # Error Correction using Coarse Grid
            for i in range(self.num_iter[1]):
                coarse = BVPSolver(self.bvp, n - 1, src_grid=self.residual().coarsen(), num_iter=self.num_iter)
                coarse.solve()
                self.sol_val += coarse.sol_grid.fine().val

            # Post-Smoothing
            for i in range(self.num_iter[2]):
                self.relax()

    def solve(self):
        self.multigrid()

    def exact_error(self):
        if self.bvp.exact_sol_func is None:
            raise Exception('exact_sol_func is not specified.')
        
        exact_error = np.array([self.bvp.exact_sol_func(xi) for xi in self.sol_grid.x]) - self.sol_val
        return Grid(self.bvp.domain, self.n, val=exact_error)