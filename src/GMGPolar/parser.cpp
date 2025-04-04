#include "../../include/GMGPolar/gmgpolar.h"

/* Specifies whether user input is required */
enum
{
    OPTIONAL = 0,
    REQUIRED = 1
};

void GMGPolar::parseGrid()
{
    R0_                 = parser_.get<double>("R0");
    Rmax_               = parser_.get<double>("Rmax");
    nr_exp_             = parser_.get<int>("nr_exp");
    ntheta_exp_         = parser_.get<int>("ntheta_exp");
    anisotropic_factor_ = parser_.get<int>("anisotropic_factor");
    divideBy2_          = parser_.get<int>("divideBy2");
    write_grid_file_    = parser_.get<int>("write_grid_file") != 0;
    load_grid_file_     = parser_.get<int>("load_grid_file") != 0;
    file_grid_radii_    = parser_.get<std::string>("file_grid_radii");
    file_grid_angles_   = parser_.get<std::string>("file_grid_angles");
    DirBC_Interior_     = parser_.get<int>("DirBC_Interior") != 0;
}

void GMGPolar::parseGeometry()
{
    const int alphaValue = parser_.get<int>("alpha_coeff");
    if (alphaValue == static_cast<int>(AlphaCoeff::POISSON) ||
        alphaValue == static_cast<int>(AlphaCoeff::SONNENDRUCKER) || alphaValue == static_cast<int>(AlphaCoeff::ZONI) ||
        alphaValue == static_cast<int>(AlphaCoeff::ZONI_SHIFTED)) {
        alpha_ = static_cast<AlphaCoeff>(alphaValue);
    }
    else {
        throw std::runtime_error("Invalid alpha coefficient value.\n");
    }

    const int problemValue = parser_.get<int>("problem");
    if (problemValue == static_cast<int>(ProblemType::CARTESIAN_R2) ||
        problemValue == static_cast<int>(ProblemType::CARTESIAN_R6) ||
        problemValue == static_cast<int>(ProblemType::POLAR_R6) ||
        problemValue == static_cast<int>(ProblemType::REFINED_RADIUS)) {
        problem_ = static_cast<ProblemType>(problemValue);
    }
    else {
        throw std::runtime_error("Invalid choice for the problem.\n");
    }

    const int geometryValue = parser_.get<int>("geometry");
    if (geometryValue == static_cast<int>(GeometryType::CIRCULAR) ||
        geometryValue == static_cast<int>(GeometryType::SHAFRANOV) ||
        geometryValue == static_cast<int>(GeometryType::CZARNY) ||
        geometryValue == static_cast<int>(GeometryType::CULHAM)) {
        geometry_ = static_cast<GeometryType>(geometryValue);
    }
    else {
        throw std::runtime_error("Invalid choice for the geometry\n");
    }

    const int betaValue = parser_.get<int>("beta_coeff");
    if (betaValue == static_cast<int>(BetaCoeff::ZERO) || betaValue == static_cast<int>(BetaCoeff::ALPHA_INVERSE)) {
        beta_ = static_cast<BetaCoeff>(betaValue);
    }
    else {
        throw std::runtime_error("Invalid beta coefficient value.\n");
    }

    alpha_jump_ = parser_.get<double>("alpha_jump");

    kappa_eps_ = parser_.get<double>("kappa_eps");
    delta_e_   = parser_.get<double>("delta_e");

    selectTestCase();
}

void GMGPolar::parseMultigrid()
{
    FMG_                     = parser_.get<int>("FMG") != 0;
    FMG_iterations_          = parser_.get<int>("FMG_iterations");
    const int FMG_cycleValue = parser_.get<int>("FMG_cycle");
    if (FMG_cycleValue == static_cast<int>(MultigridCycleType::V_CYCLE) ||
        FMG_cycleValue == static_cast<int>(MultigridCycleType::W_CYCLE) ||
        FMG_cycleValue == static_cast<int>(MultigridCycleType::F_CYCLE)) {
        FMG_cycle_ = static_cast<MultigridCycleType>(FMG_cycleValue);
    }
    else {
        throw std::runtime_error("Invalid extrapolation value.\n");
    }
    const int extrapolationValue = parser_.get<int>("extrapolation");
    if (extrapolationValue == static_cast<int>(ExtrapolationType::NONE) ||
        extrapolationValue == static_cast<int>(ExtrapolationType::IMPLICIT_EXTRAPOLATION) ||
        extrapolationValue == static_cast<int>(ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING) ||
        extrapolationValue == static_cast<int>(ExtrapolationType::COMBINED)) {
        extrapolation_ = static_cast<ExtrapolationType>(extrapolationValue);
    }
    else {
        throw std::runtime_error("Invalid extrapolation value.\n");
    }
    max_levels_           = parser_.get<int>("maxLevels");
    pre_smoothing_steps_  = parser_.get<int>("preSmoothingSteps");
    post_smoothing_steps_ = parser_.get<int>("postSmoothingSteps");

    const int cycleValue = parser_.get<int>("multigridCycle");
    if (cycleValue == static_cast<int>(MultigridCycleType::V_CYCLE) ||
        cycleValue == static_cast<int>(MultigridCycleType::W_CYCLE) ||
        cycleValue == static_cast<int>(MultigridCycleType::F_CYCLE)) {
        multigrid_cycle_ = static_cast<MultigridCycleType>(cycleValue);
    }
    else {
        throw std::runtime_error("Invalid multigrid cycle value.\n");
    }

    max_iterations_ = parser_.get<int>("maxIterations");

    const int normValue = parser_.get<int>("residualNormType");
    if (normValue == static_cast<int>(ResidualNormType::EUCLIDEAN) ||
        normValue == static_cast<int>(ResidualNormType::WEIGHTED_EUCLIDEAN) ||
        normValue == static_cast<int>(ResidualNormType::INFINITY_NORM)) {
        residual_norm_type_ = static_cast<ResidualNormType>(normValue);
    }
    else {
        throw std::runtime_error("Invalid residual norm type.\n");
    }

    double absTol = parser_.get<double>("absoluteTolerance");
    if (absTol < 0) {
        absolute_tolerance_ = std::nullopt;
    }
    else {
        absolute_tolerance_ = absTol;
    }
    double relTol = parser_.get<double>("relativeTolerance");
    if (relTol < 0) {
        relative_tolerance_ = std::nullopt;
    }
    else {
        relative_tolerance_ = relTol;
    }
}

void GMGPolar::parseGeneral()
{
    verbose_         = parser_.get<int>("verbose");
    paraview_        = parser_.get<int>("paraview") != 0;
    max_omp_threads_ = parser_.get<int>("maxOpenMPThreads");
    omp_set_num_threads(max_omp_threads_);
    thread_reduction_factor_                 = parser_.get<double>("threadReductionFactor");
    const int stencilDistributionMethodValue = parser_.get<int>("stencilDistributionMethod");
    if (stencilDistributionMethodValue == static_cast<int>(StencilDistributionMethod::CPU_TAKE) ||
        stencilDistributionMethodValue == static_cast<int>(StencilDistributionMethod::CPU_GIVE)) {
        stencil_distribution_method_ = static_cast<StencilDistributionMethod>(stencilDistributionMethodValue);
    }
    else {
        throw std::runtime_error("Invalid stencil distribution method.\n");
    }
    cache_density_profile_coefficients_ = parser_.get<int>("cacheDensityProfileCoefficients") != 0;
    cache_domain_geometry_              = parser_.get<int>("cacheDomainGeometry") != 0;
}

void GMGPolar::initializeGrid()
{
    parser_.add<double>("R0", 'r', "Interior radius of the disk.", OPTIONAL, 1e-5);
    parser_.add<double>("Rmax", 'R', "Exterior radius of the disk.", OPTIONAL, 1.3);
    parser_.add<int>("nr_exp", 'n', "Number of nodes (exponents) in the radial direction.", OPTIONAL, 5);
    parser_.add<int>("ntheta_exp", '\0', "Number of nodes (exponents) in the angular direction.", OPTIONAL, -1);
    parser_.add<int>("anisotropic_factor", '\0', "Defines anisotropic discretization in the radial direction.",
                     OPTIONAL, 0);
    parser_.add<int>("divideBy2", '\0', "Refines the grid globally `divideBy2` times.", OPTIONAL, 0);
    parser_.add<int>("write_grid_file", '\0', "Enable writing the finest PolarGrid to a file.", OPTIONAL, 0,
                     cmdline::oneof(0, 1));
    parser_.add<int>("load_grid_file", '\0', "Enable loading the finest PolarGrid from a file.", OPTIONAL, 0,
                     cmdline::oneof(0, 1));
    parser_.add<std::string>("file_grid_radii", '\0',
                             "Path to the file containing radii values for grid divisions in the radial direction.",
                             OPTIONAL, "");
    parser_.add<std::string>("file_grid_angles", '\0',
                             "Path to the file containing theta values for grid divisions in the angular direction.",
                             OPTIONAL, "");
    parser_.add<int>(
        "DirBC_Interior", '\0',
        "Defines the boundary condition on the interior circle: Across-origin (0), Dirichlet boundary (1).", OPTIONAL,
        0, cmdline::oneof(0, 1));
}

void GMGPolar::initializeGeometry()
{
    parser_.add<int>("geometry", '\0',
                     "Defines the shape of the cross-section: Circular (0), Shafranov (1), Czarny (2), Culham (3).",
                     OPTIONAL, 0, cmdline::oneof(0, 1, 2, 3));
    parser_.add<double>("alpha_jump", '\0', "Defines the radius of rapid decay of the density profile alpha.", OPTIONAL,
                        0.0);
    parser_.add<double>("kappa_eps", 'k', "Defines the elongation of the geometry.", OPTIONAL, 0.0);
    parser_.add<double>("delta_e", 'd', "Defines the outward radial displacement of the center of flux.", OPTIONAL,
                        0.0);

    parser_.add<int>("problem", '\0',
                     "Defines the problem to solve (exact solution): CartesianR2 (0), PolarR6 (1), CartesianR6 (2), "
                     "RefinedRadius (3).",
                     OPTIONAL, 0, cmdline::oneof(0, 1, 2, 3));
    parser_.add<int>("alpha_coeff", '\0',
                     "Defines the alpha coefficient: Poisson (0), Sonnendrucker (1), Zoni (2), Zoni-Shifted (3).",
                     OPTIONAL, 1, cmdline::oneof(0, 1, 2, 3));
    parser_.add<int>("beta_coeff", '\0', "Defines the beta coefficient: beta=0 (0), beta=1/alpha (1).", OPTIONAL, 0,
                     cmdline::oneof(0, 1));
}

void GMGPolar::initializeMultigrid()
{
    parser_.add<int>("FMG", '\0', "Specifies whether the initial approximation is obtained by nested iteration.",
                     OPTIONAL, 0, cmdline::oneof(0, 1));
    parser_.add<int>("FMG_iterations", '\0', "Specifies the number of FMG iterations.", OPTIONAL, 2);
    parser_.add<int>("FMG_cycle", '\0', "Specifies the type of FMG Cycle: V-cycle (0), W-cycle (1), F-cycle (2).",
                     OPTIONAL, 0, cmdline::oneof(0, 1, 2));
    parser_.add<int>("extrapolation", 'e',
                     "Specifies the type of extrapolation: No extrapolation (0), Implicit extrapolation (1), Implicit "
                     "extrapolation with full grid smoothing (2), Combination of both methods (3).",
                     OPTIONAL, 0, cmdline::oneof(0, 1, 2, 3));
    parser_.add<int>("maxLevels", 'l', "Defines the maximum number of levels in the multigrid scheme.", OPTIONAL, -1);
    parser_.add<int>("preSmoothingSteps", '\0', "Number of pre-smoothing steps.", OPTIONAL, 1);
    parser_.add<int>("postSmoothingSteps", '\0', "Number of post-smoothing steps.", OPTIONAL, 1);
    parser_.add<int>("multigridCycle", '\0', "Type of Multigrid Cycle: V-cycle (0), W-cycle (1), F-cycle (2).",
                     OPTIONAL, 0, cmdline::oneof(0, 1, 2));

    parser_.add<int>("residualNormType", '\0',
                     "Type of Residual Norm: Euclidean (0), Weighted Euclidean (1), Infinity (2).", OPTIONAL, 0,
                     cmdline::oneof(0, 1, 2));
    parser_.add<int>("maxIterations", '\0', "Maximum number of Multigrid iterations.", OPTIONAL, 150);
    parser_.add<double>("absoluteTolerance", '\0', "Convergence achieved when absolute tolerance is reached.", OPTIONAL,
                        1e-8);
    parser_.add<double>("relativeTolerance", '\0', "Convergence achieved when relative tolerance is reached.", OPTIONAL,
                        1e-8);
}

void GMGPolar::initializeGeneral()
{
    parser_.add<int>(
        "verbose", '\0',
        "Controls the verbosity of the output. Higher values produce more detailed diagnostic information.", OPTIONAL,
        1);
    parser_.add<int>("paraview", '\0', "Specifies whether to generate Paraview output files.", OPTIONAL, 0);
    parser_.add<int>("maxOpenMPThreads", '\0', "Defines the maximum number of OpenMP threads used.", OPTIONAL, 1);
    parser_.add<double>("threadReductionFactor", '\0', "Reduction factor for threads at coarser levels.", OPTIONAL,
                        1.0);
    parser_.add<int>("stencilDistributionMethod", '\0',
                     "Specifies how to distribute the stencil: CPU_Take (0), CPU_Give (1).", OPTIONAL, 0,
                     cmdline::oneof(0, 1));
    parser_.add<int>("cacheDensityProfileCoefficients", '\0',
                     "Specifies whether to cache the density profile coefficients or recalculate them dynamically.",
                     OPTIONAL, 1, cmdline::oneof(0, 1));
    parser_.add<int>("cacheDomainGeometry", '\0',
                     "Specifies whether to cache the domain geometry or recalculate it dynamically.", OPTIONAL, 1,
                     cmdline::oneof(0, 1));
}
