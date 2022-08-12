# Multigrid in Various Languages

## Newtonian Spherically Symmetric Star
Newtonian gravitational potential of spherically symmetric star is governed by
$$\frac{d^{2}\Phi}{dr^{2}}+\frac{2}{r}\frac{d\Phi}{dr}=4\pi\rho$$
where $r$ is the radius of spherical coordinate, $\Phi\left(r\right)$ is the gravitational potential, and $\rho\left(r\right)$ is the mass density.

## Density Profile and Boundary Condition
We give a density profile as
$$\rho\left(r\right)=\begin{cases} \rho_{c}\left(1-\frac{r^{2}}{r_{\mathrm{s}}^{2}}\right) &: r\leq r_{\mathrm{s}} \\\\ 0 &: r>r_{\mathrm{s}} \end{cases}$$
where $\rho_{\mathrm{c}}$ is the mass density at the center and $r_{\mathrm{s}}$ is the radius of star. Therefore, we need boundary conditions as
$$\begin{cases} \frac{d\Phi}{dr}=0 &: r=0 \\\\ \Phi=0 &: r=\infty   \end{cases}$$

## Coordinate Transformation

## Finite Difference Method
