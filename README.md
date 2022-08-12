# Multigrid in Various Languages

## Newtonian Spherically Symmetric Star
Newtonian gravitational potential of spherically symmetric star is governed by
$$\frac{d^{2}\Phi}{dr^{2}}+\frac{2}{r}\frac{d\Phi}{dr}=4\pi\rho$$
where $r$ is the radius of spherical coordinate, $\Phi\left(r\right)$ is the gravitational potential, and $\rho\left(r\right)$ is the mass density.

## Density Profile and Boundary Condition
We give a density profile as
$$\rho\left(r\right)=\begin{cases} \rho_{c}\left(1-\frac{r^{2}}{r_{\mathrm{s}}^{2}}\right) &: r< r_{\mathrm{s}} \\\\ 0 &: r\geq r_{\mathrm{s}} \end{cases}$$
where $\rho_{\mathrm{c}}$ is the mass density at the center and $r_{\mathrm{s}}$ is the radius of star. Therefore, we give boundary conditions at the center and infinity as
$$\begin{cases} \frac{d\Phi}{dr}=0 &: r=0 \\\\ \Phi=0 &: r=\infty   \end{cases}$$

## Coordinate Transformation
We transform the coordinate $r\in\left[0,\infty\right)$ to $s\in\left[0,1\right)$ by
$$s=\frac{r/r_{\mathrm{s}}}{1+r/r_{\mathrm{s}}}$$
Then, we have the equation
$$\frac{d^{2}\Phi}{ds^{2}}+\frac{2}{s}\frac{d\Phi}{ds}=4\pi\rho \frac{r_{\mathrm{s}}^{2}}{\left(1-s\right)^{4}}$$
with the density profile
$$\rho\left(s\right)=\begin{cases} \rho_{c} \left\\{ 1-\left(\frac{s}{1-s}\right)^{2} \right\\} &: 0 \leq s <\frac{1}{2} \\\\ 0 &: \frac{1}{2} \leq s < 1 \end{cases}$$
and the boundary conditions
$$\begin{cases} \frac{d\Phi}{ds}=0 &: s=0 \\\\ \Phi=0 &: s=1   \end{cases}$$


## Finite Difference Method
Let us consider a discretizaion of the domain $[0,1]$ with uniform $N$ cells. Its adjacent points have coordinates given by
$$s_{i}=ih$$
where $h=1/N$ and $i\in\left\\{0,1,\cdots,N\right\\}$. We descritize derivatives in the equation and the boundary conditions as
$$\frac{d^{2}\Phi}{ds^{2}}+\frac{2}{s}\frac{d\Phi}{ds}=\frac{\Phi_{i+1} + \Phi_{i-1} - 2\Phi_{i}}{h^{2}}+\frac{2}{ih}\frac{\Phi_{i+1}-\Phi_{i-1}}{2h}+O\left(h^{2}\right)$$
$$\frac{d\Phi}{ds}=\frac{1}{h}\left(-\frac{3}{2}\Phi_{0}+2\Phi_{1}-\frac{1}{2}\Phi_{2}\right)+O\left(h^{2}\right)$$
where $\Phi_{i}=\Phi\left(s_{i}\right)$.

## Relaxation Method
