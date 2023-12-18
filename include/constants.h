/* 
* Copyright (C) 2019-2023 The GMGPolar Development Team
*
* Authors: Philippe Leleux, Christina Schwarz, Martin J. Kühn, Carola Kruse, Ulrich Rüde
*
* Contact: 
*    Carola Kruse <Carola.Kruse@CERFACS.fr>
*    Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

/*!
 * \file constants.h
 * \brief Header defining constants and enumerations of Controls/Info
 * \author M. Kuehn, C. Kruse, P. Leleux
 * \version 0.0
 */
#ifndef CONSTANTS_HXX_
#define CONSTANTS_HXX_

#include <cmath>

/// Defines the control parameters indices in a safe way
namespace Param
{
//! To be used with the gyro::icntl vector.
enum icontrols
{
    /*! \brief (WIP) Optimized code
         *
         * 0: old version (non-optimized version for validation in the beginning, not maintained for some time; do not use)
         * 1: new version (default)
         */
    optimized,
    /*! \brief Verbose level
         *
         * Defines the level of verbose outputs. For higher levels all output from lower levels is always given.
         * 0: No output
         * 1: Minimal convergence output.
         * 2: Minimal problem setting details, iteration output and timings.
         * 3: Print all parameters and more detailed problem settings.
         * 4: Information on called functions.
         * 5: Prints information on multigrid levels.
         * 6: Print out all kind of variable values.
         */
    verbose,
    /*! \brief openmp paralellization
         *
         * Defines the number of threads we use
         */
    openmp,
    /*! \brief Matrix-free implementation
         *
         * Defines if a matrix free implementation is used
         * (A is applied) or not (A is assembled)
         */
    matrix_free,
    /*! \brief Check the implementation
         *
         * Compare the results of the current implementation with
         * structures extracted from the MATLAB implementation.
         * Only the setup is run, no multigrid cycle.
         */
    debug,
    /*! \brief Number of nodes (exponents)
         *
         * Defines the number of nodes in each direction:
         * - nr = 2^(nr_exp-1)
         * - ntheta = 2^(ceil(log2(nr))-1)
         *
         * More detailed:
         *
         * First, the number of nodes in r-direction is computed as follows: 
         * without any ansisotropy (fac_ani=0), we have nr=2^(nr_exp - 1) equally distributed nodes;
         * with anisotropy: we fisrt create nr=2^(nr_exp) - 1^(fac_ani) equally distributed nodes, and then, the grid is refined fac_ani-times;
         * Lastly, regardless of using any anisotropy or not, we refine again by splitting all intervals at the midpoint 
         * (to obtain a uniform grid refinement between the two finest levels for a successful extrapolation).
         *
         * Second, the number of nodes in theta-direction is computed as follows: 
         * ntheta= 2^(ceil(log_2(nr)))
         *
         */
    nr_exp,
    /*! \brief The parameter ntheta_exp is not used in our simulations. We generally define the number of theta intervals similar to the number of intervals in r.
         */
    ntheta_exp,
    /*! \brief Anisotropic discretization in 'edge' region for r
         *
         * Defines if we use anisotropic discretization in r-direction
         *
         * Possible values are:
         * - -1: input list of r and theta coordinates
         * - 0: no anisotropy
         * - >0: anisotropy (automatic) (number of refinements)
         *
         */
    fac_ani,
    /*! \brief Not used
         * (we never have an ansisotropy in theta-direction, so the value is always 0/false)
         */
    theta_aniso,
    /*! \brief Smoothing steps
         *
         * Number of pre- (v1) and post- (v2) -smoothing steps
         */
    v1,
    v2,
    /*! \brief Type of MG cycle
         *
         * Type of multigrid cycle:
         * - 1: V-cycle (default setting)
         * - 2: W-cycle
         */
    cycle,
    /*! \brief Defines the form of the considered cross-section: Circular or stretched geometry. If `mod_pk=0`, we consider a circular geometry. If `mod_pk>1`, it always goes with a particular choice of `kappa_eps` and `delta_e` to describe realistic Tokamak cross sections. 
    
    For more details, we refer to:
    - Bouzat, N., Bressan, C., Grandgirard, V., Latu, G., Mehrenberger, M.: Targeting Realistic Geometry in Tokamak Code Gysela. (2018)
    - Zoni, E., Güçlü, Y.: Solving hyperbolic-elliptic problems on singular mapped disk-like
domains with the method of characteristics and spline finite elements.  (2019)
    - Bourne et al.: Solver comparison for Poisson-like equations on tokamak geometries. (2023)
         *
         * Defines the shape of the geometry:
         * - 0: kappa_eps = delta_e = 0 (circular geometry)
         * - 1: kappa_eps=0.3, delta_e=0.2 (stretched and deformed circle, also denoted Shafranov geometry)
         * - 2: kappa_eps=0.3, delta_e=1.4 (stretched and deformed circle, also denoted Czarny geometry)
         */
    mod_pk,
    /*! \brief Compute rho
         *
         * Defines if we compute rho, i.e. the reduction factor for the
         * residual through the iteration.
         */
    compute_rho,
    /*! \brief Maximum level for MG
         *
         * Defines the maximum level used for the multigrid scheme:
         * -1: limit such that there are 3 points in r, and 4 points in theta on the finest grid.
         */
    level,
    /*! \brief Plot or not
         *
         * Plot or not
         */
    plotit,
    /*! \brief Compare with exact solution
         *
         * Compute the the exact solution without MG for comparison
         * (ok with manufactured solutions).
         */
    solveit,
    /*! \brief Maximum iterations for MG
         *
         * Maximum number of iterations for the multigrid scheme.
         */
    maxiter,
    /*! \brief Periodic boundary condition in theta
         *
         * Periodic boundary condition in theta. The other possibility was
         * Dirichlet BC but we do not consider this anymore.
         */
    periodic,
    /*! \brief Specifies if theta has a periodicity condition
         *
         * Specifies if theta has a periodicity condition
         */
    origin_NOT_coarse,
    /*! \brief Smoother (3 and 13 only)
         *
         * Defines the smoother we use:
         * - 0: PointRB (not implemented)
         * - 1: CircZGS (not implemented)
         * - 2: RadZGS (not implemented)
         * - 3: AltZGS (coupled circle-radial version) (default setting)
         * - 13: Decoupled AltZGS (C-R: BJ, W-B: GS) (faster to execute but leads to more iterations)
         * - 4: optimizedAltZGS (not implemented)
         * - 5: MixZGS (not implemented)
         */
    smoother,
    /*! \brief Discretization scheme (3 only)
         *
         * Defines the discretization scheme we use
         * 0: FD (useless)
         * 1: 5p (not implemented)
         * 2: 7p (not implemented)
         * 3: 9p
         * 4: P1 (not implemented)
         * 5: P1nonst (not implemented)
         * 6: P2 (not implemented)
         */
    discr,
    /*! \brief Extrapolation
         *
         * 0: no extrapolation
         * 1: implicit extrapolation with adapted smoothing on finest grid (default setting)
         * 2: experimental version of implicit extrapolation with full grid smoothing (residual stopping criterion not functional, use with care)
         * 3: extrap_integer (FE only)
         * 4: no_extrap_integr (FE only)
         */
    extrapolation,
    /*! \brief Boundary conditions on the interior circle
         *
         * - 0: Across the origin
         * - 1: Dirichlet
         */
    DirBC_Interior,
    /*! \brief Output paraview file
         *
         * Output paraview file for visualization
         */
    paraview,
    /*! \brief Divide the intervals of the grid by 2^(divideBy2)
         *
         * Divides the intervals of the grid by 2^(divideBy2)
         * (defines how often to split the intervals of the grid at the midpoint)
         */
    divideBy2,
    /*! \brief Problem to solve
         *
         * Defines the problem to solve
         * - 5: the classical problem (Zoni2019, aligned with the cartesian grid)
         * - 6: more realistic problem (Emily & Virginie, aligned with the polar grid)
         */
    prob,
    /*! \brief Coefficient alpha
         *
         * Defines the coefficient alpha
         * 0: arctan coeff (Sonnendrucker2019, private communication)
         * 1: tanh coeff (Zoni2019)
         * 2: tanh coeff (Zoni2019 with steeper and shifted jump)
         */
    alpha_coeff,
    /*! \brief Coefficient beta
         *
         * Defines the coefficient
         * 0: beta = 0
         * 1: beta: 1/alpha for some cases, different in others (see coeffs in test_cases)
         */
    beta_coeff,
    /*! \brief Norm for stopping criterion
         * Defines the norm used for the residual in the stopping criterion
         * 0: L2 norm scaled by initial res.
         * 1: Infinitiy norm scaled by initial res.
         * 2: L2 norm
         * 3: Infinitiy norm
         */
    res_norm,
    /*! \brief Read or write radii/angles files
         * When f_grid_r/f_grid_theta not empty,
         * 0: read radii/angles file for the grid.
         * 1: write radii/angles file after creating the grid.
         */
    write_radii_angles,
    /*! \brief Check the error
         * Compute the error w.r.t. the theoretical solution
         * 0: do note check
         * 1: check. If the parameter "sol_in" is given, the error is computed w.r.t. a solution read from a file.
         */
    check_error,
};
//! To be used with the gyro::dcntl vector.
enum dcontrols
{
    /*! \brief Radius of the disk
         *
         * Interior (R0) and exterior (R) radius of the disk-like shape
         */
    r0_DB,
    R0,
    R,
    /*! \brief angles inside the disk
         *
         * Start/end angle for the disk-like shape discretization
         */
    THETA0,
    THETA,
    /*! \brief Parameters of the grid
         *
         * Defines the shape of the geometry:
         * - if mod_pk = 0: circular geometry
         *       [x=r.cos(theta), y=r.sin(theta)]
         * - elif mod_pk=1: stretched circular geometry (Shafranov geometry)
         *       [x=(1-kappa_eps)r*cos(theta)-delta_e*r^2, y=(1+kappa_eps)r*sin(theta)]
         * - elif mod_pk=2: stretched circular geometry (Czarny geometry)
         *       [x=1/kappa_eps ( 1 - sqrt{1 + kappa_eps ( kappa_eps + 2r*cos(theta) )} ), y=y_0 + (delta_e \xi r sin(theta)) / (1 + kappa_eps x(r, theta))]
         * [For more details, we refer to Bourne et al. (2023).]
         * with
         *     - kappa_eps: Elongation (For the Shafranov geometry, the parameter is denoted kappa, for the Czarny geometry, it is epsilon.
         *     - delta_e: Shafranov shift (outward radial
         *       displacement of the centre of flux) (For the Shafranov geometry, this parameter is denoted delta, for the Czarny geometry, it is e.)
         */
    kappa_eps,
    delta_e,
    /*! \brief Distance to boundary
         *
         * Defines the distance below which nodes are considered part of the boundary.
         */
    tol_bound_check,
    /*! \brief Threshold for convergence
         *
         * Threshold for convergence with stopping criterion:
         * - scaled residual for extrapolation < 2
         * - evolution of error w.r.t. theoretical solution for extrapolation 2
         */
    rel_red_conv,
    /*! \brief Timings
         */
    t_arr_art_att, // Evaluation of arr, art, and att
    t_coeff, // Evaluation of alpha and beta (subcounter of t_arr_art_att)    
    t_sol,
    t_detDFinv,
    t_trafo,
};
//! To be used with the gyro::info vector.
enum info
{
};
//! To be used with the gyro::dinfo vector.
enum dinfo
{
};
//! Stencils are represented by vectors with:
// - -1: not in the stencil
// - n: position of the entry w.r.t other entries of the same node
enum stencil
{
    bottom_left,
    bottom,
    bottom_right,
    left,
    middle,
    right,
    top_left,
    top,
    top_right,
};
} // namespace Param

//const double PI = 3.141592653589793238463;
const double PI = M_PI;

// See Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
// used as Param::mod_pk ("modified polar coordinates") in GMGPolar
enum geometry_type
{
    CIRCULAR   = 0, // simple circular domain
    SHAFRANOV  = 1, // Fig. 6a
    TRIANGULAR = 2, // Fig. 6b (also denoted Czarny)
    CULHAM     = 3 // Fig. 18
};

enum alpha_val
{
    SONNENDRUCKER = 0,
    ZONI          = 1,
    ZONI_SHIFTED  = 2,
    POISSON       = 3,
};

// Defines the manufactured solution to compare the computed error against.
// see Kuehn et al. https://doi.org/10.1007/s10915-022-01802-1
// or Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
enum problem_type
{
    FLAT           = 1,
    REFINED_RADIUS = 4,
    CARTESIAN_R2   = 5, // 
    POLAR_R6       = 6, // Bourne et al., Eq. (22)
    CARTESIAN_R6   = 7, // Bourne et al., Eq. (23)
};
#endif // CONSTANTS_HXX
