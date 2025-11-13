

## Ordinary Viscosity

Alvelius Forcing: Determines energy injection rate $P$ and the
forcing scale $k_f$.

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

## Hyper Viscosity

The incompressible Navier-Stokes equation with hyperviscosity is given by

$$ \partial_t \mathbf{u} + \mathbf{u} \cdot \nabla \mathbf{u} = - \nabla p + (-1)^{{h+1} } \nu \Delta^h \mathbf{u} + \mathbf{f} $$

For the case with ordinary viscosity we set h = 1.

The dissipation scale for hyperviscous turbulence is given by

$$\eta = \left( \frac{\nu^3}{\epsilon} \right) ^{ 1 / 6h - 2 }$$

For $h = 1$ we arrive at the already knwon Kolmogorov scale.

Generalizing the above considerations, the viscosity in the hyperviscous case must be set as

$$ \boxed{ \nu = \left( \frac{c_{res}}{k_{max} } P^{1/m} \right)^{m/3} }$$

, with $m = 6h - 2$.