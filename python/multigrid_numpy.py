import numpy as np

from multigrid import BVPSolver, Grid
from smart_slice import SmartSlice, cast


class BVPSolver_numpy(BVPSolver):
    """BVPSolver which uses SmartSlice and its associated array. It is meant to
    reuse the same function that are used with BVPSolver, but with faster
    execution.
    """

    def relax(self):
        # Aliases
        h = self.h
        x = self.x
        N = self.N
        sol_val = self.sol_val
        src_val = self.src_val

        # Red Sweep
        sol_val[0] = self.bvp.relax_left_func(sol_val, src_val, x[0], h)
        i = SmartSlice(2, N-1, 2)
        cast(sol_val)[i] = self.bvp.relax_middle_func(
            cast(sol_val), cast(src_val), cast(x)[i], h, i)
        sol_val[N] = self.bvp.relax_right_func(sol_val, src_val, x[N], h)

        # Black Sweep
        i = SmartSlice(1, N, 2)
        cast(sol_val)[i] = self.bvp.relax_middle_func(
            cast(sol_val), cast(src_val), cast(x)[i], h, i)

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
        i = SmartSlice(1, N, 1)
        cast(res)[i] = self.bvp.res_middle_func(
            cast(sol_val), cast(src_val), cast(x)[i], h, i)

        # Residual at the rightmost point
        res[N] = self.bvp.res_left_func(sol_val, src_val, x[N], h)

        return Grid(self.bvp.domain, self.n, val=res)
