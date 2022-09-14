import numpy as np
from multigrid import BVPSolver, Grid
from numba_relax import relax_redblack_inplace, residual_jacobian


class BVPSolver_numba(BVPSolver):
    """BVP solver using numba-version of relax and residual method. While it takes
    BVP instance for initialization, `domain`, `src_func` and `exact_sol_func`
    atriutes are only used.

    """

    def relax(self):
        # Aliases
        h = self.h
        x = self.x
        sol_val = self.sol_val
        src_val = self.src_val

        relax_redblack_inplace(src_val, sol_val, x, h)

    def residual(self):
        # Aliases
        h = self.h
        x = self.x
        sol_val = self.sol_val
        src_val = self.src_val

        res = residual_jacobian(src_val, sol_val, x, h)

        return Grid(self.bvp.domain, self.n, val=res)
