from scipy.interpolate import interp1d
from copy import deepcopy
import numpy as np

def is_positive_integer(x):
    return x > 0 and isinstance(x, int)

def is_non_negative_integer(x):
    return x >= 0 and isinstance(x, int)

def get_n(x):
    return int(np.log2(x - 1))

def is_p2_1(x):
    '''
        check if x is 2**n + 1 where n is non-negative integer.
    '''
    return isinstance(x, int) and x > 1 and x == 2**get_n(x) + 1

class Grid:
    def super__init__(self, x1, x2, val=None, shape=None, size=None):
        '''
            x1: float
            x2: float
            val: numpy.ndarray of float with shape (..., size)
            shape: (..., size)
        '''
        if not x1 < x2:
            raise Exception('x1 should be smaller than x2.')
        self.x1 = x1
        self.x2 = x2
        if val is None:
            if shape is None:
                if size is None:
                    raise Exception('size is required when val and shape are None')
                else:
                    shape = (size,)
            try:
                np.ndarray(shape)
            except:
                raise Exception('shape should be tuple of non-negative integers.')
            val = np.zeros(shape)

        self.size = val.shape[-1]
        self.set_val(val)
        self.x = np.linspace(x1, x2, num=self.size)
        self.h = (x2 - x1)/self.size

    def set_val(self, val):
        '''
            val: numpy.ndarray of float with shape (..., size)
        '''
        val = np.array(val)
        if val.shape[-1] != self.size:
            raise Exception('shape of val should be (..., {}).'.format(self.size))
        self.val = val

    def rms(self):
        return np.linalg.norm(self.val, ord=2)/np.sqrt(self.size)

    def __init__(self, x1, x2, val=None, shape=None, n=None):
        if val is None:
            if shape is None or shape[-1] is None:
                if n is None:
                    raise Exception('n is required when shape has None.')
                shape = () if shape is None else shape
                shape = shape[:-1] + (2**n + 1,)
            else:
                if not is_p2_1(shape[-1]):
                    raise Exception('shape should be (..., 2**n + 1).')
        else:
            if not is_p2_1(val.shape[-1]):
                raise Exception('shape of val should be (..., 2**n + 1).')

        self.super__init__(x1, x2, val=val, shape=shape)
        self.n = get_n(self.size)

    def fine(self, kind='linear'):
        finer = Grid(self.x1, self.x2, shape=self.val.shape[:-1]+(None,), n=self.n+1)
        finer.set_val(interp1d(self.x, self.val, kind=kind)(finer.x))
        return finer

    def coarsen(self):
        val = self.val
        coarser = Grid(self.x1, self.x2, shape=self.val.shape[:-1]+(None,), n=self.n-1)
        coarser.val[..., 0] = val[..., 0]
        coarser.val[..., -1] = val[..., -1]
        for i in range(1, coarser.size - 1):
            coarser.val[..., i] = 1/4*val[..., 2*i - 1] + 1/2*val[..., 2*i] + 1/4*val[..., 2*i + 1]
        return coarser

class BVP:
    def __init__(self, x1, x2, relax_left_func, relax_middle_func, relax_right_func, res_left_func, res_middle_func, res_right_func, src_func, exact_sol_func=None):
        self.x1 = x1
        self.x2 = x2
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
        self.bvp = bvp
        self.n = n
        self.sol_grid = Grid(bvp.x1, bvp.x2, n=n)
        self.src_grid = Grid(bvp.x1, bvp.x2, val=np.array([bvp.src_func(xi) for xi in self.sol_grid.x])) if src_grid is None else src_grid
        self.num_iter = num_iter

        self.h = self.sol_grid.h
        self.x = self.sol_grid.x
        self.N = self.sol_grid.size - 1
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

        return Grid(self.sol_grid.x1, self.sol_grid.x2, val=res)

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
        return Grid(self.sol_grid.x1, self.sol_grid.x2, val=exact_error)