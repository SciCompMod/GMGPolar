#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::parseGrid() {
    R0 = parser_.get<scalar_t>("R0");
    Rmax = parser_.get<scalar_t>("R");
    nr_exp = parser_.get<int>("nr_exp");
    ntheta_exp = parser_.get<int>("ntheta_exp");
    anisotropic_factor = parser_.get<int>("fac_ani");
    divideBy2 = parser_.get<int>("divideBy2");
    write_grid_file = parser_.get<int>("write_grid_file") != 0;
    load_grid_file = parser_.get<int>("load_grid_file") != 0;
    file_grid_r = parser_.get<std::string>("file_grid_r");
    file_grid_theta = parser_.get<std::string>("file_grid_theta");
}

void GMGPolar::parseGeometry() {
    DirBC_Interior = parser_.get<int>("DirBC_Interior") != 0;

    const int alphaValue = parser_.get<int>("alpha_coeff");
    if (alphaValue == SONNENDRUCKER || alphaValue == ZONI || alphaValue == ZONI_SHIFTED || alphaValue == POISSON)
        alpha = static_cast<alpha_coeff>(alphaValue);
    else throw std::runtime_error("Invalid alpha coefficient value.\n");

    const int betaValue = parser_.get<int>("beta_coeff");
    if (betaValue == ZERO || betaValue == ALPHA_INVERSE)
        beta = static_cast<beta_coeff>(betaValue);
    else throw std::runtime_error("Invalid beta coefficient value.\n");

    const int problemValue = parser_.get<int>("problem");
    if (problemValue == CARTESIAN_R2 || problemValue == POLAR_R6 || problemValue == CARTESIAN_R6 || problemValue == REFINED_RADIUS)
        problem = static_cast<problem_type>(problemValue);
    else throw std::runtime_error("Invalid choice for the problem.\n");

    const int geometryValue = parser_.get<int>("geometry");
    if (geometryValue == CIRCULAR || geometryValue == SHAFRANOV || geometryValue == TRIANGULAR || geometryValue == CULHAM)
        geometry = static_cast<geometry_type>(geometryValue);
    else throw std::runtime_error("Invalid choice for the geometry\n");

    kappa_eps = parser_.get<scalar_t>("kappa_eps");
    delta_e = parser_.get<scalar_t>("delta_e");
}

void GMGPolar::parseMultigrid() {
    extrapolation = parser_.get<int>("extrapolation");
    maxLevels = parser_.get<int>("maxLevels");
    v1 = parser_.get<int>("v1");
    v2 = parser_.get<int>("v2");
    cycle = parser_.get<int>("cycle");
}

void GMGPolar::parseGeneral() {
    maxOpenMPThreads = parser_.get<int>("maxOpenMPThreads");
    finestLevelThreads = parser_.get<int>("finestLevelThreads");
    threadReductionFactor = parser_.get<scalar_t>("threadReductionFactor");
    omp_set_num_threads(maxOpenMPThreads);
    verbose = parser_.get<int>("verbose");
}

void GMGPolar::initializeGrid() {
    // Disk Shape Parameters
    // Defines the interior (R0) and exterior (R) radii of the disk-like shape.
    parser_.add<scalar_t>(
        "R0", 'r', "Interior radius of the disk", 
        OPTIONAL, 1e-5
    );
    parser_.add<scalar_t>(
        "R", 'R', "Exterior radius of the disk", 
        OPTIONAL, 1.3
    );
    // Defines the number of nodes in each direction:
    // - nr = 2^(nr_exp-1)
    // - ntheta = 2^(ceil(log2(nr))-1)
    //
    // More detailed:
    // First, the number of nodes in r-direction is computed as follows: 
    // without any anisotropy (fac_ani=0), we have nr=2^(nr_exp - 1) equally distributed nodes;
    // with anisotropy: we first create nr=2^(nr_exp) - 1^(fac_ani) equally distributed nodes, and then, the grid is refined fac_ani-times;
    // Lastly, regardless of using any anisotropy or not, we refine again by splitting all intervals at the midpoint 
    // (to obtain a uniform grid refinement between the two finest levels for a successful extrapolation).
    //
    // Second, the number of nodes in theta-direction is computed as follows: 
    // ntheta= 2^(ceil(log_2(nr)))
    parser_.add<int>(
        "nr_exp", 'n', 
        "Number of nodes (exponents) in the radial direction.", 
        OPTIONAL, 4
    );
    // Note: The parameter ntheta_exp is not used in our simulations.
    // We generally define the number of theta intervals similar to the number of intervals in r.
    parser_.add<int>(
        "ntheta_exp", '\0', 
        "Number of nodes (exponents) in the angular direction.", 
        OPTIONAL, -1
    );
    // Anisotropic discretization in the 'edge' region for r.
    // Defines whether anisotropic discretization is applied in the r-direction.
    // Possible values:
    // -  0: No anisotropy
    // - >0: Anisotropy (automatic) (number of refinements)
    parser_.add<int>(
        "fac_ani", '\0', 
        "Defines anisotropic discretization in r-direction.", 
        OPTIONAL, 0
    );
    // Divide the intervals of the grid by 2^(divideBy2).
    // Defines how often to split the intervals of the grid at the midpoint.
    parser_.add<int>(
        "divideBy2", '\0', 
        "Refines the grid globally divideBy2 times.", 
        OPTIONAL, 0
    );
    // Toggle whether to write the finest PolarGrid to a file
    parser_.add<int>(
        "write_grid_file", '\0', "Enable writing the finest PolarGrid to a file.",
        OPTIONAL, 0, cmdline::oneof(0,1)
    );
    // Toggle whether to load the finest PolarGrid from a file
    parser_.add<int>(
        "load_grid_file", '\0', "Enable loading the finest PolarGrid from a file.",
        OPTIONAL, 0, cmdline::oneof(0,1)
    );
    // Loading the input finest PolarGrid from a list of r and theta coordinates
    parser_.add<std::string>(
        "file_grid_r", '\0', "Path to the file containing radii values for grid divisions in the r-direction.",
        OPTIONAL, ""
    );
    parser_.add<std::string>(
        "file_grid_theta", '\0', "Path to the file containing theta values for grid divisions in the theta-direction.",
        OPTIONAL, ""
    );
}

void GMGPolar::initializeGeometry() {

    parser_.add<int>(
        "DirBC_Interior", '\0', "Defines the boundary condition on the interior circle. Across-origin(0), Dirichlet-boundary(1).",
        OPTIONAL, 0, cmdline::oneof(0,1)
    );

    // Defines the coefficient alpha
    // 0: arctan coeff (Sonnendrucker2019, private communication)
    // 1: tanh coeff (Zoni2019)
    // 2: tanh coeff (Zoni2019 with steeper and shifted jump)
    // 3: ??? coeff (Poisson)
    parser_.add<int>(
        "alpha_coeff", '\0', "Defines the coefficient alpha. Sonnendrucker(0), Zoni(1), Zoni-Shifted(2), Poisson(3).",
        OPTIONAL, 0, cmdline::oneof(0,1,2,3)
    );
    // Defines the coefficient beta
    // 0: beta = 0
    // 1: beta: 1/alpha for some cases, different in others (see coeffs in test_cases)
    parser_.add<int>(
        "beta_coeff", '\0', "Defines the coefficient beta. beta=0 (0), beta=1/alpha (1).",
        OPTIONAL, 0, cmdline::oneof(0, 1)
    );
    // Defines the problem to solve
    // 5: the classical problem (Zoni2019, aligned with the cartesian grid)
    // 6: more realistic problem (Emily & Virginie, aligned with the polar grid)
    parser_.add<int>(
        "problem", '\0', "Defines the problem to solve. CartesianR2 (0), PolarR6 (1), CartesianR6 (2), RefinedRadius (3).",
        OPTIONAL, 0, cmdline::oneof(0,1,2,3)
    );
    // Defines the form of the considered cross-section: 
    // Circular or stretched geometry. If `mod_pk=0`, we consider a circular geometry. 
    // If mod_pk > 1, it always goes with a particular choice of `kappa_eps` and `delta_e` 
    // to describe realistic Tokamak cross sections. 
    //
    //For more details, we refer to:
    // - Bouzat, N., Bressan, C., Grandgirard, V., Latu, G., Mehrenberger, M.: 
    //   Targeting Realistic Geometry in Tokamak Code Gysela. (2018)
    // - Zoni, E., Güçlü, Y.:
    //   Solving hyperbolic-elliptic problems on singular mapped disk-like domains 
    //   with the method of characteristics and spline finite elements.  (2019)
    // - Bourne et al.: 
    //   Solver comparison for Poisson-like equations on tokamak geometries. (2023)
    // 
    // Defines the shape of the geometry:
    // - 0: kappa_eps = delta_e = 0 (circular geometry)
    // - 1: kappa_eps=0.3, delta_e=0.2 (stretched and deformed circle, also denoted Shafranov geometry)
    // - 2: kappa_eps=0.3, delta_e=1.4 (stretched and deformed circle, also denoted Czarny geometry)
    parser_.add<int>( // CHANGED FROM MOD_PK TO GEOMETRY
        "geometry", '\0', "Defines the form of the considered cross-section. Circular (0), Shafranov (1), Czarny/Triangular (2), Culham (3)",
        OPTIONAL, 0, cmdline::oneof(0,1,2,3)
    );
    // Defines the shape of the geometry:
    // - if mod_pk = 0: circular geometry
    //      [x=r.cos(theta),
    //       y=r.sin(theta)]
    // - elif mod_pk=1: stretched circular geometry (Shafranov geometry)
    //      [x=(1-kappa_eps)r*cos(theta)-delta_e*r^2, 
    //       y=(1+kappa_eps)r*sin(theta)]
    // - elif mod_pk=2: stretched circular geometry (Czarny geometry)
    //      [x=1/kappa_eps ( 1 - sqrt{1 + kappa_eps ( kappa_eps + 2r*cos(theta) )} ),
    //       y=y_0 + (delta_e \xi r sin(theta)) / (1 + kappa_eps x(r, theta))]
    // - MISSING: MOD_PK=3 EXPLAINATION
    //
    // [For more details, we refer to Bourne et al. (2023).]
    //
    // - kappa_eps: Elongation (For the Shafranov geometry), denoted by kappa
    //              Elongation (For the Czarny geometry), denoted by epsilon
    // 
    // - delta_e: Shafranov shift (outward radial displacement of the centre of flux) 
    //            For the Shafranov geometry, this parameter is denoted delta.
    //            For the Czarny geometry, it is denoted by e.
    parser_.add<scalar_t>(
        "kappa_eps", 'k', "Defines the Elongation of the geometry.",
        OPTIONAL, 0.0
    );
    parser_.add<scalar_t>(
        "delta_e", 'd', "Defines the outward radial displacement of the centre of flux.",
        OPTIONAL, 0.0
    );
}

void GMGPolar::initializeMultigrid() {
    // Extrapolation options:
    //  * 0: No extrapolation
    //  * 1: Implicit extrapolation with adapted smoothing on the finest grid (default setting)
    //  * 2: Experimental version of implicit extrapolation with full grid smoothing (residual stopping criterion not functional, use with care)
    //  * 3: Extrapolation using integers (Finite Element only)
    //  * 4: No extrapolation using integers (Finite Element only)
    parser_.add<int>(
        "extrapolation", 'e', 
        "Specifies the type of extrapolation method to be used.", 
        OPTIONAL, 0, cmdline::oneof(0,1,2,3,4)
    );
    // Maximum levels for multigrid scheme:
    //  * -1: Limit such that there are 3 points in radial direction and 4 points in theta direction on the finest grid.
    parser_.add<int>(
        "maxLevels", 'l', 
        "Defines the maximum number of levels used in the multigrid scheme.", 
        OPTIONAL, -1
    );
    // Number of Pre-Smoothing steps
    parser_.add<int>(
        "v1", '\0', 
        "Defines the number of pre-smoothing steps.", 
        OPTIONAL, 1
    );
    // Number of Post-Smoothing steps
    parser_.add<int>(
        "v2", '\0', 
        "Defines the number of post-smoothing steps.", 
        OPTIONAL, 1
    );
    // Type of Multigrid Cycle
    // - 1: V-cycle (default setting)
    // - 2: W-cycle
    parser_.add<int>(
        "cycle", '\0', 
        "Type of Multigrid Cycle.", 
        OPTIONAL, 1
    );
}

void GMGPolar::initializeGeneral() {
    // OpenMP paralellization
    parser_.add<int>(
        "maxOpenMPThreads", '\0', 
        "Defines the maximum number of OpenMP threads used.", 
        REQUIRED, 1
    );
    parser_.add<int>(
        "finestLevelThreads", '\0', 
        "Optimal number of OpenMP threads for the finest level, can exceed maxOpenMPThreads.", 
        REQUIRED, 1
    );
    parser_.add<scalar_t>(
        "threadReductionFactor", '\0', 
        "Thread reduction factor to coarser levels.", 
        OPTIONAL, 1.0
    );

    // Defines the level of verbose outputs.
    // For higher levels all output from lower levels is always given.
    // 0: No output
    // 1: Minimal convergence output.
    // 2: Minimal problem setting details, iteration output and timings.
    // 3: Print all parameters and more detailed problem settings.
    // 4: Information on called functions.
    // 5: Prints information on multigrid levels.
    // 6: Print out all kind of variable values.
    parser_.add<int>(
        "verbose", '\0', 
        "Defines the level of verbose outputs.", 
        OPTIONAL, 2, cmdline::oneof(0,1,2,3,4,5,6)
    );
}
