# SpecTurb

A pseudo-spectral code for 3D forced MHD turbulence simulations on parallel CPU architectures.

## What does it do?

*SpecTurb* solves the incompressible MHD equations

$$  
\begin{align*}
\partial_t \mathbf{u} &= - \mathbf{u} \cdot \nabla \mathbf{u} - \nabla p + \nu \Delta \mathbf{u} + ( \nabla \times \mathbf{B} ) \times \mathbf{B} + \mathbf{f}, \\  
\partial_t \mathbf{B} &= \nabla \times ( \mathbf{u} \times \mathbf{B} ) + \eta \Delta \mathbf{B}, \\
\nabla \cdot \mathbf{u} &= 0, \\
\nabla \cdot \mathbf{B} &= 0.
\end{align*}
$$

on a periodic domain.

To drive the system into an equilibrium state, we apply the forcing $\mathbf{f}$ proposed by [Alvelius](https://doi.org/10.1063/1.870050).  
With this forcing, we can set the energy production rate $P$ as well as the forced wavenumber $k_f$.

The macroscopic scales, where the energy is injected into the fluid, are defined as

- $T$: Large-eddy turnover time  
- $L$: Integral length scale  
- $U$: Large-eddy velocity  

and are related by

$$ U = \frac{L}{T}. $$

Since the energy injection rate is connected to the macroscopic scales by

$$ P = \frac{U^2}{T} = \frac{L^2}{T^3}, $$

we can determine $T$ and $L = 2 \pi / k_f$ from the forcing and adjust the injection rate $P$ accordingly.  
To obtain a convenient timescale for the simulations, we set $T = 1$.

The energy dissipation scale $\eta$ is given by 

$$ \eta = \frac{\nu^{3/4}}{\epsilon^{1/4}}, $$

where $\nu$ denotes the (kinematic) viscosity of the fluid and $\epsilon$ is the energy dissipation rate.

For fully developed turbulence, the energy input equals the output:

$$ \epsilon = P, $$

and the dissipation scale $\eta$ is determined by the resolution $N$ of the simulation.

The maximum wavenumber $k_{\max}$, when applying the 2/3 de-aliasing rule, is given by

$$ k_{\max} = \frac{2}{3} \frac{N}{2} = \frac{N}{3}, $$

if the total domain size is set to $2\pi$.

A common resolution requirement is of the form

$$ \eta \, k_{\max} = c_{\text{res}}, $$

where a good resolution is said to be achieved for $c_{\text{res}} = 3$, and a sufficient resolution for $c_{\text{res}} = 1.5$.  
We usually apply the latter to maximize the energy cascade.

The resolution requirement thus determines the viscosity as

$$ \boxed{ \nu = \left( \frac{c_{\text{res}}}{k_{\max}} P^{1/4} \right)^{4/3} }. $$

To obtain a wider energy spectrum, it is possible to also simulate hyperviscous fluids, i.e. we solve

$$  
\begin{align*}
\partial_t \mathbf{u} &= - \mathbf{u} \cdot \nabla \mathbf{u} - \nabla p + (-1)^{h+1}\nu \Delta^h \mathbf{u} + ( \nabla \times \mathbf{B} ) \times \mathbf{B} + \mathbf{f}, \\  
\partial_t \mathbf{B} &= \nabla \times ( \mathbf{u} \times \mathbf{B} ) + (-1)^{h+1}\eta \Delta^h \mathbf{B}, \\
\nabla \cdot \mathbf{u} &= 0, \\
\nabla \cdot \mathbf{B} &= 0.
\end{align*}
$$

For the case of ordinary viscosity, we set $h = 1$.

The dissipation scale for hyperviscous turbulence is given by

$$ \eta = \left( \frac{\nu^3}{\epsilon} \right)^{1 / (6h - 2)}. $$

For $h = 1$, we recover the classical Kolmogorov scale.

Generalizing the above considerations, the viscosity in the hyperviscous case must be set as

$$ \boxed{ \nu = \left( \frac{c_{\text{res}}}{k_{\max}} P^{1/m} \right)^{m/3} }, $$

with $m = 6h - 2$.

---

## How to use it?

1. Clone the repository using either HTTPS or SSH:

```bash
# HTTPS
git clone https://github.com/MikeWilbert/specTurb.git

# SSH
git clone git@github.com:MikeWilbert/specTurb.git
```

2. Specify your parameters in the config file 'CONFIG/config.json'

3. Make sure your MPI compiler and the FFTW3 library are correctly linked in the Makefile and run `make`

4. Build and run the code with:

`mpiexec -np [#processes] SpecTurb`

## How to setup?

The following parameters can be specified in the config file 'CONFIG/config.json':

| Group      | Parameter  | Description                                                                                                                                                                                                      |
| ---------- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| resolution | N          | Spatial resolution in each of the three directions.                                                                                                                                                              |
|            | px, py     | Resolution of the 2D processor grid. **Note:** `(px × py)` must match the `-np` parameter when running with MPI.                                                                                                 |
| turbulence | Prm        | Magnetic Prandtl number $Pr_m = \nu/\eta$.                                                                                                                                                                       |
|            | hyp        | Hyperviscosity order $h$.                                                                                                                                                                                        |
|            | c_res      | Resolution of the micro scales $c_{\text{res}}$.                                                                                                                                                                 |
| setup      | which      | Specifies one of the following setups: (0) zero initial conditions, (1) Orszag–Tang vortex, (2) small turbulent fields (white noise with power-law slope), (3) read binary file specified by `setup:binary_dir`. |
|            | binary_dir | Directory to read a binary file containing initial fields (used with option `setup:which = 3`).                                                                                                                  |
| forcing    | on         | Boolean value enabling or disabling the forcing.                                                                                                                                                                 |
|            | k          | Mean forced mode $k_f$.                                                                                                                                                                                          |
|            | dk         | Forcing interval, i.e. the force is applied to modes in $[k - dk, k + dk]$.                                                                                                                                      |
| B_bg       | on         | Boolean value enabling a constant background magnetic field.                                                                                                                                                     |
|            | dE         | Expected magnetic energy of the turbulent field.                                                                                                                                                                 |
|            | E0_dE      | Energy ratio of the background magnetic field to the turbulent field.                                                                                                                                            |
| output     | dir        | Path to the output directory.                                                                                                                                                                                    |
|            | interval   | Output interval measured in large-eddy turnover times.                                                                                                                                                           |
|            | end        | End time of the simulation (in large-eddy turnover times).                                                                                                                                                       |
| restart    | dir        | Restart directory. If not specified (i.e. `dir = ""`), the output directory is used.                                                                                                                             |
|            | step       | Output step to restart from. If set to 0, no restart is performed and the specified setup is used.                                                                                                               |


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
