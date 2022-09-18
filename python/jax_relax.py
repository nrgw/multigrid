from jax.config import config
config.update("jax_enable_x64", True)
from jax import jit
import jax.numpy as jnp


def relax_redblack_inplace(src, sol, x, h):

    sol = jnp.array(sol)
    # Red Sweep
    sol = sol.at[0].set(sol[1] - src[0]*h**2/6)
    s = x[2:-1:2]
    v = (sol[3::2]*(1 + h/s)
         + sol[1:-2:2]*(1 - h/s)
         - src[2:-1:2]*h**2
         )/2
    sol = sol.at[2:-1:2].set(v)
    sol = sol.at[-1].set(0)

    # Black Sweep
    s = x[1:-1:2]
    v = (sol[2::2]*(1 + h/s)
         + sol[:-2:2]*(1 - h/s)
         - src[1:-1:2]*h**2
         )/2
    sol = sol.at[1:-1:2].set(v)

    return sol

def relax_jacobian(src, sol, x, h):

    relaxed = jnp.zeros_like(src)

    relaxed = relaxed.at[0].set(sol[1] - src[0]*h**2/6)
    s = x[1:-1]
    v = (sol[2:]*(1 + h/s)
         + sol[:-2]*(1 - h/s)
         - src[1:-1]*h**2
         )/2
    relaxed = relaxed.at[1:-1].set(v)
    relaxed = relaxed.at[-1].set(0)

    return relaxed

def residual_jacobian(src, sol, x, h):
    # Aliases

    res = jnp.zeros_like(src)
    relaxed = relax_jacobian(src, sol, x, h)

    # Residual at the leftmost point
    res_left = -(relaxed[0]-sol[0])*6/h**2
    res = res.at[0].set(res_left)
    # Residual at middle points
    res_middle = -(relaxed[1:-1]-sol[1:-1])*2/h**2
    res = res.at[1:-1].set(res_middle)

    # Residual at the rightmost point
    res = res.at[-1].set(0)

    return res

jrelax_redblack_inplace = jit(relax_redblack_inplace)
jresidual_jacobian = jit(residual_jacobian)
