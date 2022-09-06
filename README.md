#  Multigrid Codes in Various Languages

We implement the multigrid method in various languages to compare performance. The languages we used are:
* Fortran
* C
* C++
* Python
* Julia
* Mathematica
* Rust

## Problem Setting

### Spherically Symmetric Star

The gravitational potential is governed by the Poisson's equation:

$$
\nabla^2 \Phi=4\pi G \rho,
$$

where $\Phi$, $\rho$ and $G$ are the gravitational potential, mass density and gravitational constant, respectively.
Therefore, the gravitational potential due to a spherically symmetric object can be written in the spherical coordinate:

$$
\frac{d^{2}\Phi}{dr^{2}}+\frac{2}{r}\frac{d\Phi}{dr}=4\pi \rho,
$$

where $r$ is the radius of spherical coordinate, $\Phi\left(r\right)$ is the gravitational potential, $\rho\left(r\right)$ is the mass density, and we set $G=1$.

### Density Profile and Boundary Condition

We assume the a density profile profile follows a parabolic function in radius as follows:

$$
    \rho\left(r\right)=\begin{cases}
        \rho_{\mathrm{c}}\left(1-\frac{r^{2}}{r_{\mathrm{s}}^{2}}\right) &: r< r_{\mathrm{s}}
        \\\\ 0 &: r\geq r_{\mathrm{s}}
    \end{cases},
$$

where $\rho_{\mathrm{c}}$ is the mass density at the center and $r_{\mathrm{s}}$ is the radius of star. We can easily impose two boundary conditions at the center (using symmetry) and infinity (asymptotic behavior of gravitational potential) as

$$
    \begin{cases}
        \frac{d\Phi}{dr}=0 &: r=0
        \\\\ \Phi=0 &: r=\infty
    \end{cases}.
$$

### Coordinate Transformation
To cover the spatial infinity ( $r=\infty$ ) in the computational domain, several special coordinate transformation functions are ommonly use for the coordinate compactification.
We use the following special function:

$$
    s=\frac{r/r_{\mathrm{s}}}{1+r/r_{\mathrm{s}}},
$$

where it transform from $r\in\left[0,\infty\right)$ to $s\in\left[0,1\right)$.
Then, we can recast the Poisson's equation in spherical coordinate:

$$
    \frac{d^{2}\Phi}{ds^{2}}+\frac{2}{s}\frac{d\Phi}{ds}=4\pi\rho \frac{r_{\mathrm{s}}^{2}}{\left(1-s\right)^{4}},
$$

with the boundary conditions

$$
    \begin{cases}
        \frac{d\Phi}{ds}=0 &: s=0
        \\\\ \Phi=0 &: s=1
    \end{cases}.
$$

And the density profile becomes

$$
    \rho\left(r\right)=\begin{cases}
        \rho_{\mathrm{c}}\left[1-\left(\frac{s}{1-s}\right)^2\right] &: r< r_{\mathrm{s}}
        \\\\ 0 &: r\geq r_{\mathrm{s}}
    \end{cases}.
$$

### Finite Difference Method
Let us consider a discretizaion of the domain $[0,1]$ with uniform $N$ cells. Their adjacent points have coordinates given by
$$s_{i}=ih,$$
where $h=1/N$ and $i\in\left\\{0,1,\cdots,N\right\\}$. The subscript $i$ represents $i$-th grid point starting from $0$. We descritize LHS of the Poisson's equation as follows:

$$
    \frac{d^{2}\Phi}{ds^{2}}+\frac{2}{s}\frac{d\Phi}{ds}=\frac{\Phi_{i+1} + \Phi_{i-1} - 2\Phi_{i}}{h^{2}}+\frac{2}{ih}\frac{\Phi_{i+1}-\Phi_{i-1}}{2h}+O\left(h^{2}\right)
$$

And at $s=0$, the above equation can be written as
$$3\frac{d^2\Phi}{ds^2}=6\frac{\Phi_1-\Phi_0}{h^2}+O\left(h^{2}\right).$$
Here, we use the L'Hospital's Rule to regularize the term $\frac{2}{s}\frac{d\Phi}{ds}$ at $s=0$ since $\frac{d\Phi}{ds}=0$ at $s=0$, and the reflective boundary condition ( $\Phi_{-1}=\Phi_{1}$ ) is also applied.

### Relaxation Method
For $i=N$, we fix

$$
    \Phi_{N}=0.
$$

For $i=0$,

$$
    \Phi_{0}\leftarrow\Phi_{1}-\frac{1}{6}\left(4\pi\rho_{\mathrm{c}}r_{\mathrm{s}}^{2}\right) h^2.
$$

For $i\in\left\\{1,2,\cdots,N-1\right\\}$,

$$
    \Phi_{i}\leftarrow\frac{1}{2}\left[\left(1+\frac{1}{i}\right)\Phi_{i+1}+\left(1-\frac{1}{i}\right)\Phi_{i-1}-\left(4\pi\rho_{i}\frac{r_{\mathrm{s}}^{2}}{\left(1-s_{i}\right)^{4}}\right)h^{2}\right].
$$

## Multigrid Method
