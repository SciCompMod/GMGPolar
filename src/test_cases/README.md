The following files define the geometry, coefficients of the ordinary differential equation and the manufactured solution. We always solve

$-\nabla \cdot \left( \alpha \nabla u \right) + \beta u = f\quad \text{ in }\quad \Omega$

The naming is defined in `constants.h` and it follows the following rules:

- The first part, e.g., `PolarR6` or `CartesianR6`, defines the manufactured solution.
- The **potential** second part `Gyro` or empty, defines if the $\beta$ coefficient is nontrivial or zero.
- The following part, e.g., `Sonnendrucker` or `Poisson`, defines the value of the $\alpha$ coefficient.
- The last part, e.g., `Shafranov` or `Triangular`, defines the geometry to be used

For the particular solutions, see:
- Bourne et al. - Solver comparison for Poisson-like equations on tokamak geometries (2023) https://doi.org/10.1016/j.jcp.2023.112249 
- Kuehn, Kruse, Ruede - Implicitly extrapolated geometric multigrid on disk-like domains for the gyrokinetic Poisson equation from fusion plasma applications (2022) https://doi.org/10.1007/s10915-022-01802-1

A great thanks to Emily Bourne for preparing these scripts!