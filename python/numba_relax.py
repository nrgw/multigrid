import numpy as np
from numba import stencil, njit

@stencil(neighborhood = ((-1, 1),))
def relax_middle(src, sol, x, h):
    return (
        + sol[1]*(1 + h/x[0])
        + sol[-1]*(1 - h/x[0])
        - src[0]*h**2
    )/2

# This is for red-black. We assume that src has slice of [i:j:2] while the sol
# has a slice of [i-1:j:2]. Therefore, original indexing of sol[i-1] becomes
# sol[0] and sol[i+1] becoms sol[1]
@stencil(neighborhood = ((0, 0),))
def relax_middle2(src2, sol2m1, x2, h):
    return (
        + sol2m1[1]*(1 + h/x2[0])
        + sol2m1[0]*(1 - h/x2[0])
        - src2[0]*h**2
    )/2

@njit
def relax_jacobian(src, sol, x, h):

    relaxed = relax_middle(src, sol, x, h)
    relaxed[0] = sol[1] - src[0]*x[0]**2/6
    relaxed[-1] = 0

    return relaxed

@njit
def relax_redblack_inplace(src, sol, x, h):
    # sol = sol.copy()

    sol[0] = sol[1] - src[0]*x[0]**2/6

    relax_middle2(src[2::2], sol[1::2], x[2::2], h, out=sol[2::2])

    sol[-1] = 0.

    relax_middle2(src[1::2], sol[0::2], x[1::2], h, out=sol[1::2])

    return sol

@njit
def residual_jacobian(src, sol, x, h):
    res = np.zeros_like(src)

    relaxed = relax_jacobian(src, sol, x, h)

    res[0] = -(relaxed[0] - sol[0])*6/h**2
    res[1:-1] = -(relaxed[1:-1] - sol[1:-1])*2/h**2
    res[-1] = 0.

    return res
