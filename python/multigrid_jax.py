import numpy as np
from multigrid import BVPSolver, Grid
from jax_relax import jrelax_redblack_inplace, jresidual_jacobian


class BVPSolver_jax(BVPSolver):
    """BVP solver using jax-version of relax and residual method. While it takes
    BVP instance for initialization, `domain`, `src_func` and `exact_sol_func`
    atriutes are only used.
    """

    def relax(self):
        # Aliases
        h = self.h
        x = self.x
        sol_val = self.sol_val
        src_val = self.src_val

        v = jrelax_redblack_inplace(src_val, sol_val, x, h)
        self.sol_val[:] = v

    def residual(self):
        # Aliases
        h = self.h
        x = self.x
        sol_val = self.sol_val
        src_val = self.src_val

        # using the returned jnp array from jresidual somehow fails.
        res = np.zeros(sol_val.shape)
        _res = jresidual_jacobian(src_val, sol_val, x, h)
        res[:] = _res

        return Grid(self.bvp.domain, self.n, val=res)
