from scipy.interpolate import interp1d
import numpy as np
from smart_slice import SmartSlice, cast

class Grid:
    SHAREDMEM = {}
    DOMAINX = {}

    @classmethod
    def preallocate_memory(cls, n, domain):
        x1, x2 = domain
        for m in range(1, n+1):
            cls.SHAREDMEM[(2**m+1,)] = [np.zeros((2**m+1,)) for i in range(3)]
            cls.DOMAINX[(domain, 2**m+1)] = np.linspace(x1, x2, 2**m+1)

    def __init__(self, domain, n, val=None, func=None, shape=None):
        self.domain = domain
        self.n = n

        N = 2**n
        x1 = domain[0]
        x2 = domain[1]

        self.N = N
        self.x = self.DOMAINX[(domain, N+1)]

        self.h = (x2 - x1)/N
        self._memory_borrowed = False
        
        if val is not None:
            if val.shape[0] != N + 1:
                raise Exception('Shape of val is invalid.')
            self.val = val
        elif func is not None:
            self.val = np.array([func(xi) for xi in self.x])
        elif shape is not None:
            self.val = self.SHAREDMEM[(self.N + 1,) + shape].pop()
            self.val.fill(0)
            self._memory_borrowed = True
        else:
            raise Exception('One of val, func, and shape should be specified.')
            
    def __del__(self):
        if self._memory_borrowed:
            self.SHAREDMEM[self.val.shape].append(self.val)

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

        # coarser.val[1:coarser.N] = (1/4*val[1:2*coarser.N - 1:2] +
        #                             1/2*val[2:2*coarser.N:2] +
        #                             1/4*val[3:2*coarser.N + 1:2])

        i = SmartSlice(2, 2*coarser.N, 2)
        val = cast(val)
        coarser.val[1:coarser.N] = (1/4*val[i-1] +
                                    1/2*val[i] +
                                    1/4*val[i+1])

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

        # Aliases
        self.h = self.sol_grid.h
        self.x = self.sol_grid.x
        self.N = self.sol_grid.N
        self.sol_val = self.sol_grid.val
        self.src_val = self.src_grid.val

    def relax(self):
        # Aliases
        h = self.h
        x = self.x
        N = self.N
        sol_val = self.sol_val
        src_val = self.src_val

        # Red Sweep
        sol_val[0] = self.bvp.relax_left_func(sol_val, src_val, x[0], h)
        for i in range(2, N, 2):
            sol_val[i] = self.bvp.relax_middle_func(sol_val, src_val, x[i], h, i)
        sol_val[N] = self.bvp.relax_right_func(sol_val, src_val, x[N], h)

        # Black Sweep
        for i in range(1, N, 2):
            sol_val[i] = self.bvp.relax_middle_func(sol_val, src_val, x[i], h, i)

    def residual(self):
        # Aliases
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
        # Pre-Smoothing
        for i in range(self.num_iter[0]):
            self.relax()
        if self.n != 1:
            # Error Correction using Coarse Grid
            for i in range(self.num_iter[1]):
                coarse = type(self)(self.bvp, self.n - 1, src_grid=self.residual().coarsen(), num_iter=self.num_iter)
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
