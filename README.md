# SpecTurb

A pseudo-spectral code producing forced MHD turbulence in 3D.

## What does it do?

*SpecTurb* solves the incompressible MHD equations

$$  
\begin{align*}
\partial_t \mathbf{u} &= - \mathbf{u} \cdot \nabla \mathbf{u} - \nabla p + \nu \Delta \mathbf{u} + ( \nabla \times \mathbf{B} ) \times \mathbf{B} + \mathbf{f} \\  
\partial_t \mathbf{B} &= \nabla \times ( \mathbf{u} \times \mathbf{B} ) + \eta \Delta \mathbf{B}\\
\nabla \cdot \mathbf{u} &= 0 \\
\nabla \cdot \mathbf{B} &= 0
\end{align*}
$$

on a periodic domain.

To drive the system into an equilibrium state, we apply the forcing $\mathbf{f}$ propoesed by [Alvelius](https://doi.org/10.1063/1.870050). With this forcing, we can set the energy production rate $P$ as well as the forced wave length $k_f$.

The macroscopic scales, where the energy is injected into the fluid are called

- $T$ : Large-Eddy Turnover Time
- $L$ : Inertial Length Scale
- $U$ : Large Eddy Velocity

, where the 3 scales are connected by

$$ U = L/T $$

Since the energy injection rate is connted to the macroscopic scales by

$$P = \frac{U^2}{T} = \frac{L^2}{T^3}$$

we can simply determine $T$ and $L = 2 \pi / k_f$ by the forcing and adjust the injection rate $P$ accordingly. To get a reasonable timescale for the simulations, we set $T = 1$.

The energy dissipation scale $\eta$ is given by 

$$\eta = \frac{\nu^{3/4}}{\epsilon^{1/4}}$$

Here, $\nu$ denotes the (kinetic) viscosity of the fluid and $\epsilon$ is the energy dissipation rate.

For in fully developed turbulence the energy input equals the output, we have

$$\epsilon = P$$

and the dissipation scale $\eta$ is determined by the resolution $N$ of the simulation.

The maximum wavenumber $k_{max}$, when applying the 2/3 dealiasing rule, is given by

$$ k_{max} = \frac{2}{3} \frac{N}{2} = \frac{N}{3} $$

if the total domain size is set to $2 \pi$.

A common resolution requirement is often used of the form

$$ \eta \, k_{max} = c_{res} $$

, where a good resolution is said to be achieved for $c_{ref} = 3$ and a sufficient resolution is given by $c_{res} = 1.5$. We usually apply the second case to maximize the energy cascade.

The resolution requirement thus determines the viscosity to

$$ \boxed{ \nu = \left( \frac{c_{res}}{k_{max} } P^{1/4} \right)^{4/3} }$$

To get a wider energy spectrum, it is possible to also simulate hyperviscous fluids, i.e. there we solve

$$  
\begin{align*}
\partial_t \mathbf{u} &= - \mathbf{u} \cdot \nabla \mathbf{u} - \nabla p + (-1)^{h+1}\nu \Delta^h \mathbf{u} + ( \nabla \times \mathbf{B} ) \times \mathbf{B} + \mathbf{f} \\  
\partial_t \mathbf{B} &= \nabla \times ( \mathbf{u} \times \mathbf{B} ) + (-1)^{h+1}\eta \Delta^h \mathbf{B}\\
\nabla \cdot \mathbf{u} &= 0 \\
\nabla \cdot \mathbf{B} &= 0
\end{align*}
$$

For the case with ordinary viscosity we set h = 1.

The dissipation scale for hyperviscous turbulence is given by

$$\eta = \left( \frac{\nu^3}{\epsilon} \right) ^{ 1 / 6h - 2 }$$

For $h = 1$ we arrive at the already knwon Kolmogorov scale.

Generalizing the above considerations, the viscosity in the hyperviscous case must be set as

$$ \boxed{ \nu = \left( \frac{c_{res}}{k_{max} } P^{1/m} \right)^{m/3} }$$

, with $m = 6h - 2$.

## How to use it?

1. Clone the repository using either HTTPS or SSH:

```bash
# HTTPS
git clone https://github.com/MikeWilbert/specTurb.git

# SSH
git clone git@github.com:MikeWilbert/specTurb.git
```

2. Specify your parameters in the config file 'CONFIG/config.json'

3. Make sure your MPI compiler and the FFTW3 library are correctly linked in the Makefile.

4. Build and run the code with:

`mpiexec -np [#processes] SpecTurb`

### Config File


## Third-party libraries

- **OpenMPI** (BSD-style license, external)
  - Used for MPI parallelization
  - Users must install OpenMPI on their system
  - Tested with OpenMPI 4.x
  - License: https://www.open-mpi.org/software/ompi/license/

- **FFTW3** (GPL, external)
  - Used for FFTs
  - Users must install FFTW3
  - License: http://www.fftw.org/

- **nlohmann/json** (MIT License, included)
  - Header-only JSON parser included in `include/nlohmann/json.hpp`