from scipy.interpolate import interp1d
from copy import deepcopy
import numpy as np

def is_positive_integer(x):
    return x > 0 and isinstance(x, int)

def is_non_negative_integer(x):
    return x >= 0 and isinstance(x, int)

class Grid:
    def __init__(self, x1, x2, val=None, shape=None, size=None):
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

    # def copy(self, val=None):
    #     clone = deepcopy(self)
    #     if val is not None:
    #         clone.set_val(val)
    #     return clone

def get_n(x):
    return int(np.log2(x - 1))

def is_p2_1(x):
    '''
        check if x is 2**n + 1 where n is non-negative integer.
    '''
    return isinstance(x, int) and x > 1 and x == 2**get_n(x) + 1

class MultiGrid(Grid):
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

        super().__init__(x1, x2, val=val, shape=shape)
        self.n = get_n(self.size)

    def fine(self, kind='linear'):
        finer = MultiGrid(self.x1, self.x2, shape=self.val.shape[:-1]+(None,), n=self.n+1)
        finer.set_val(interp1d(self.x, self.val, kind=kind)(finer.x))
        return finer

    def coarsen(self):
        val = self.val
        coarser = MultiGrid(self.x1, self.x2, shape=self.val.shape[:-1]+(None,), n=self.n-1)
        coarser.val[..., 0] = val[..., 0]
        coarser.val[..., -1] = val[..., -1]
        for i in range(1, coarser.size - 1):
            coarser.val[..., i] = 1/4*val[..., 2*i - 1] + 1/2*val[..., 2*i] + 1/4*val[..., 2*i + 1]
        return coarser

class BVP:
    def __init__(self, sol, src, relax, residual):
        pass
            
def multigrid(sol, src, relax, residual):
    relax(sol, src)
    if sol.n != 1:
        err_coarse = MultiGrid(sol.x1, sol.x2, n=sol.n-1)
        multigrid(err_coarse, residual(sol, src).coarsen(), relax, residual)
        sol.val += err_coarse.fine().val
        relax(sol, src)

def fullmultigrid(sol, src, relax, residual):
    # for i in range(1, sol.n):
    pass
        