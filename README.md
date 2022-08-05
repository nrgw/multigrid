# Multigrid in Various Languages

## Boundary Value Problem
For a function $f\left(x\right)$, we give a boundary value problem with a equation,
$$\frac{d^{2}f}{dx^{2}}=g\left(x\right)$$
where $g\left(x\right)$ is a source term and boundary condition,
$$\begin{align} f\left(0\right) &=0 \\\\ f\left(1\right) &=0 \end{align}$$

## Lane-Emden Equation
Newtonian self-gravtating, spherically symmetric star with polytropic fluid is governed by
$$\frac{1}{\xi^{2}}\frac{d}{d\xi}\left(\xi^{2}\frac{d\theta}{d\xi}\right)+\theta^{n}=0$$
$$\begin{align} \rho&=\rho_{\mathrm{c}}\theta^{n} \\\\ P&=K\rho^{1+1/n} \\ \end{align}$$
where $\xi$ is a dimensionless radius, $\theta$ is a variable for matter, $\rho$ is the density, $\rho_{\mathrm{c}}$ is its value at center, $n$ is the polytropic index, $P$ is the pressure, and $K$ is a constant.

## Finite Difference Method
