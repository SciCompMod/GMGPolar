#include "../include/ConfigParser/config_parser.h"

ConfigParser::ConfigParser()
{
    // Initialize command-line options for general parameters
    parser_.add<int>("verbose", '\0', "Verbosity level.", OPTIONAL, 1);
    parser_.add<int>("paraview", '\0', "Generate ParaView output (0/1).", OPTIONAL, 0);
    parser_.add<int>("maxOpenMPThreads", '\0', "Max OpenMP threads.", OPTIONAL, 1);
    parser_.add<double>("threadReductionFactor", '\0', "Thread reduction factor.", OPTIONAL, 1.0);
    parser_.add<int>("DirBC_Interior", '\0', "Interior BC type (0=Across-origin, 1=Dirichlet).", OPTIONAL, 0,
                     cmdline::oneof(0, 1));
    parser_.add<int>("stencilDistributionMethod", '\0', "Stencil distribution (0=CPU_Take,1=CPU_Give)", OPTIONAL, 0,
                     cmdline::oneof(0, 1));
    parser_.add<int>("cacheDensityProfileCoefficients", '\0', "Cache density coefficients (0/1).", OPTIONAL, 1,
                     cmdline::oneof(0, 1));
    parser_.add<int>("cacheDomainGeometry", '\0', "Cache domain geometry (0/1).", OPTIONAL, 1, cmdline::oneof(0, 1));

    // Initialize command-line options for grid parameters
    parser_.add<double>("R0", 'r', "Interior radius of the disk.", OPTIONAL, 1e-5);
    parser_.add<double>("Rmax", 'R', "Exterior radius of the disk.", OPTIONAL, 1.3);
    parser_.add<int>("nr_exp", 'n', "Number of nodes (exponents) in radial direction.", OPTIONAL, 5);
    parser_.add<int>("ntheta_exp", '\0', "Number of nodes (exponents) in angular direction.", OPTIONAL, -1);
    parser_.add<int>("anisotropic_factor", '\0', "Anisotropic discretization factor in radial direction.", OPTIONAL, 0);
    parser_.add<int>("divideBy2", '\0', "Global refinement steps.", OPTIONAL, 0);

    // Initialize command-line options for multigrid parameters
    parser_.add<int>("extrapolation", 'e', "Extrapolation method (0=None,1=Implicit,2=FullGrid,3=Combined)", OPTIONAL,
                     0, cmdline::oneof(0, 1, 2, 3));
    parser_.add<int>("maxLevels", 'l', "Max multigrid levels.", OPTIONAL, -1);
    parser_.add<int>("preSmoothingSteps", '\0', "Pre-smoothing steps.", OPTIONAL, 1);
    parser_.add<int>("postSmoothingSteps", '\0', "Post-smoothing steps.", OPTIONAL, 1);
    parser_.add<int>("multigridCycle", '\0', "Cycle type (0=V,1=W,2=F).", OPTIONAL, 0, cmdline::oneof(0, 1, 2));
    parser_.add<int>("FMG", '\0', "Use Full Multigrid (0/1).", OPTIONAL, 0, cmdline::oneof(0, 1));
    parser_.add<int>("FMG_iterations", '\0', "FMG iterations.", OPTIONAL, 2);
    parser_.add<int>("FMG_cycle", '\0', "FMG cycle type (0=V,1=W,2=F).", OPTIONAL, 0, cmdline::oneof(0, 1, 2));
    parser_.add<int>("maxIterations", '\0', "Max solver iterations.", OPTIONAL, 150);
    parser_.add<int>("residualNormType", '\0', "Residual norm (0=Euclidean,1=Weighted,2=Infinity)", OPTIONAL, 0,
                     cmdline::oneof(0, 1, 2));
    parser_.add<double>("absoluteTolerance", '\0', "Absolute tolerance.", OPTIONAL, 1e-8);
    parser_.add<double>("relativeTolerance", '\0', "Relative tolerance.", OPTIONAL, 1e-8);

    // Initialize command-line options for geometry parameters
    parser_.add<int>("geometry", '\0', "Geometry type (0=Circular,1=Shafranov,2=Czarny,3=Culham)", OPTIONAL, 0,
                     cmdline::oneof(0, 1, 2, 3));
    parser_.add<double>("alpha_jump", '\0', "Radius for density decay.", OPTIONAL, 0.0);
    parser_.add<double>("kappa_eps", 'k', "Elongation parameter.", OPTIONAL, 0.0);
    parser_.add<double>("delta_e", 'd', "Radial displacement.", OPTIONAL, 0.0);
    parser_.add<int>("problem", '\0', "Problem type (0=CartesianR2,1=PolarR6,2=CartesianR6,3=RefinedRadius)", OPTIONAL,
                     0, cmdline::oneof(0, 1, 2, 3));
    parser_.add<int>("alpha_coeff", '\0', "Alpha coefficient (0=Poisson,1=Sonnendrucker,2=Zoni,3=Zoni-Shifted)",
                     OPTIONAL, 1, cmdline::oneof(0, 1, 2, 3));
    parser_.add<int>("beta_coeff", '\0', "Beta coefficient (0=Zero,1=1/alpha)", OPTIONAL, 0, cmdline::oneof(0, 1));

    parse(0, nullptr);
}

// Parses command-line arguments and initializes all parameters
bool ConfigParser::parse(int argc, char* argv[])
{

    if (argc != 0) {
        try {
            parser_.parse_check(argc, argv);
        }
        catch (const cmdline::cmdline_error& parse_error) {
            std::cerr << "Error: " << parse_error.what() << std::endl;
            std::cerr << "Usage: " << parser_.usage() << std::endl;
        }
    }

    // Parse general parameters from command-line arguments
    verbose_                 = parser_.get<int>("verbose");
    paraview_                = parser_.get<int>("paraview") != 0;
    max_omp_threads_         = parser_.get<int>("maxOpenMPThreads");
    thread_reduction_factor_ = parser_.get<double>("threadReductionFactor");
    DirBC_Interior_          = parser_.get<int>("DirBC_Interior") != 0;
    const int methodValue    = parser_.get<int>("stencilDistributionMethod");
    if (methodValue == static_cast<int>(StencilDistributionMethod::CPU_TAKE) ||
        methodValue == static_cast<int>(StencilDistributionMethod::CPU_GIVE)) {
        stencil_distribution_method_ = static_cast<StencilDistributionMethod>(methodValue);
    }
    else {
        throw std::runtime_error("Invalid stencil distribution method.");
    }
    cache_density_profile_coefficients_ = parser_.get<int>("cacheDensityProfileCoefficients") != 0;
    cache_domain_geometry_              = parser_.get<int>("cacheDomainGeometry") != 0;

    // Parse grid parameters from command-line arguments
    double R0              = parser_.get<double>("R0");
    double Rmax            = parser_.get<double>("Rmax");
    double nr_exp          = parser_.get<int>("nr_exp");
    int ntheta_exp         = parser_.get<int>("ntheta_exp");
    int anisotropic_factor = parser_.get<int>("anisotropic_factor");
    int divideBy2          = parser_.get<int>("divideBy2");

    // Parse multigrid parameters from command-line arguments
    FMG_                     = parser_.get<int>("FMG") != 0;
    FMG_iterations_          = parser_.get<int>("FMG_iterations");
    const int FMG_cycleValue = parser_.get<int>("FMG_cycle");
    if (FMG_cycleValue == static_cast<int>(MultigridCycleType::V_CYCLE) ||
        FMG_cycleValue == static_cast<int>(MultigridCycleType::W_CYCLE) ||
        FMG_cycleValue == static_cast<int>(MultigridCycleType::F_CYCLE)) {
        FMG_cycle_ = static_cast<MultigridCycleType>(FMG_cycleValue);
    }
    else {
        throw std::runtime_error("Invalid FMG cycle type.");
    }
    const int extrapolationValue = parser_.get<int>("extrapolation");
    if (extrapolationValue == static_cast<int>(ExtrapolationType::NONE) ||
        extrapolationValue == static_cast<int>(ExtrapolationType::IMPLICIT_EXTRAPOLATION) ||
        extrapolationValue == static_cast<int>(ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING) ||
        extrapolationValue == static_cast<int>(ExtrapolationType::COMBINED)) {
        extrapolation_ = static_cast<ExtrapolationType>(extrapolationValue);
    }
    else {
        throw std::runtime_error("Invalid extrapolation type.");
    }
    max_levels_           = parser_.get<int>("maxLevels");
    pre_smoothing_steps_  = parser_.get<int>("preSmoothingSteps");
    post_smoothing_steps_ = parser_.get<int>("postSmoothingSteps");
    const int cycleValue  = parser_.get<int>("multigridCycle");
    if (cycleValue == static_cast<int>(MultigridCycleType::V_CYCLE) ||
        cycleValue == static_cast<int>(MultigridCycleType::W_CYCLE) ||
        cycleValue == static_cast<int>(MultigridCycleType::F_CYCLE)) {
        multigrid_cycle_ = static_cast<MultigridCycleType>(cycleValue);
    }
    else {
        throw std::runtime_error("Invalid multigrid cycle type.");
    }
    max_iterations_     = parser_.get<int>("maxIterations");
    const int normValue = parser_.get<int>("residualNormType");
    if (normValue == static_cast<int>(ResidualNormType::EUCLIDEAN) ||
        normValue == static_cast<int>(ResidualNormType::WEIGHTED_EUCLIDEAN) ||
        normValue == static_cast<int>(ResidualNormType::INFINITY_NORM)) {
        residual_norm_type_ = static_cast<ResidualNormType>(normValue);
    }
    else {
        throw std::runtime_error("Invalid residual norm type.");
    }
    double absTol       = parser_.get<double>("absoluteTolerance");
    absolute_tolerance_ = (absTol >= 0) ? std::optional<double>(absTol) : std::nullopt;
    double relTol       = parser_.get<double>("relativeTolerance");
    relative_tolerance_ = (relTol >= 0) ? std::optional<double>(relTol) : std::nullopt;

    const int geometry_value = parser_.get<int>("geometry");
    GeometryType geometry_type;
    if (geometry_value == static_cast<int>(GeometryType::CIRCULAR) ||
        geometry_value == static_cast<int>(GeometryType::CULHAM) ||
        geometry_value == static_cast<int>(GeometryType::CZARNY) ||
        geometry_value == static_cast<int>(GeometryType::SHAFRANOV)) {
        geometry_type = static_cast<GeometryType>(geometry_value);
    }
    else {
        throw std::runtime_error("Invalid geometry type.");
    }
    double alpha_jump       = parser_.get<double>("alpha_jump");
    double kappa_eps        = parser_.get<double>("kappa_eps");
    double delta_e          = parser_.get<double>("delta_e");
    const int problem_value = parser_.get<int>("problem");
    ProblemType problem_type;
    if (problem_value == static_cast<int>(ProblemType::CARTESIAN_R2) ||
        problem_value == static_cast<int>(ProblemType::CARTESIAN_R6) ||
        problem_value == static_cast<int>(ProblemType::POLAR_R6) ||
        problem_value == static_cast<int>(ProblemType::REFINED_RADIUS)) {
        problem_type = static_cast<ProblemType>(problem_value);
    }
    else {
        throw std::runtime_error("Invalid problem type.");
    }
    AlphaCoeff alpha_type;
    const int alpha_value = parser_.get<int>("alpha_coeff");
    if (alpha_value == static_cast<int>(AlphaCoeff::POISSON) ||
        alpha_value == static_cast<int>(AlphaCoeff::SONNENDRUCKER) ||
        alpha_value == static_cast<int>(AlphaCoeff::ZONI) ||
        alpha_value <= static_cast<int>(AlphaCoeff::ZONI_SHIFTED)) {
        alpha_type = static_cast<AlphaCoeff>(alpha_value);
    }
    else {
        throw std::runtime_error("Invalid alpha coefficient.");
    }
    BetaCoeff beta_type;
    const int beta_value = parser_.get<int>("beta_coeff");
    if (beta_value == static_cast<int>(BetaCoeff::ZERO) || beta_value == static_cast<int>(BetaCoeff::ALPHA_INVERSE)) {
        beta_type = static_cast<BetaCoeff>(beta_value);
    }
    else {
        throw std::runtime_error("Invalid beta coefficient.");
    }

    // Construct PolarGrid
    double refinement_radius               = alpha_jump;
    std::optional<double> splitting_radius = std::nullopt;
    grid_ = PolarGrid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

    selectTestCase(geometry_type, problem_type, alpha_type, beta_type, Rmax, kappa_eps, delta_e, alpha_jump);

    if (verbose_ && argc != 0) {
        std::cout << "------------------------------\n";
        std::cout << "------ Problem Settings ------\n";
        std::cout << "------------------------------\n";

        if (typeid(*domain_geometry_) == typeid(CircularGeometry)) {
            std::cout << "Circular (Domain geometry)\n";
        }
        else if (typeid(*domain_geometry_) == typeid(ShafranovGeometry)) {
            std::cout << "Shafranov (Domain geometry)\n";
        }
        else if (typeid(*domain_geometry_) == typeid(CzarnyGeometry)) {
            std::cout << "Czarny (Domain geometry)\n";
        }
        else if (typeid(*domain_geometry_) == typeid(CulhamGeometry)) {
            std::cout << "Culham (Domain geometry)\n";
        }
        else {
            std::cout << "Unknown domain geometry\n";
        }

        if (exact_solution_ != nullptr) {
            // std::cout << "Exact Solution: ";
            if (typeid(*exact_solution_) == typeid(CartesianR2_CircularGeometry) ||
                typeid(*exact_solution_) == typeid(CartesianR2_CzarnyGeometry) ||
                typeid(*exact_solution_) == typeid(CartesianR2_ShafranovGeometry)) {
                std::cout << "CartesianR2 (Exact solution)\n";
            }
            else if (typeid(*exact_solution_) == typeid(CartesianR6_CircularGeometry) ||
                     typeid(*exact_solution_) == typeid(CartesianR6_CzarnyGeometry) ||
                     typeid(*exact_solution_) == typeid(CartesianR6_ShafranovGeometry)) {
                std::cout << "CartesianR6 (Exact solution)\n";
            }
            else if (typeid(*exact_solution_) == typeid(PolarR6_CircularGeometry) ||
                     typeid(*exact_solution_) == typeid(PolarR6_CulhamGeometry) ||
                     typeid(*exact_solution_) == typeid(PolarR6_CzarnyGeometry) ||
                     typeid(*exact_solution_) == typeid(PolarR6_ShafranovGeometry)) {
                std::cout << "PolarR6 (Exact solution)\n";
            }
            else if (typeid(*exact_solution_) == typeid(Refined_CircularGeometry) ||
                     typeid(*exact_solution_) == typeid(Refined_CulhamGeometry) ||
                     typeid(*exact_solution_) == typeid(Refined_CzarnyGeometry) ||
                     typeid(*exact_solution_) == typeid(Refined_ShafranovGeometry)) {
                std::cout << "Multi-Scale (Exact solution)\n";
            }
            else {
                std::cout << "Unknown exact solution\n";
            }
        }

        if (typeid(*density_profile_coefficients_) == typeid(PoissonCoefficients)) {
            std::cout << "α = 1, β = 0 (Poisson)\n";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(SonnendruckerCoefficients)) {
            std::cout << "α = Sonnendrücker, β = 0 (Profile coefficients)\n";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(SonnendruckerGyroCoefficients)) {
            std::cout << "α = Sonnendrücker, β = 1/α (Profile coefficients)\n";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(ZoniCoefficients)) {
            std::cout << "α = Zoni, β = 0 (Profile coefficients)\n";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(ZoniGyroCoefficients)) {
            std::cout << "α = Zoni, β = 1/α  (Profile coefficients)\n";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(ZoniShiftedCoefficients)) {
            std::cout << "α = Zoni-Shifted, β = 0 (Profile coefficients)\n";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(ZoniShiftedGyroCoefficients)) {
            std::cout << "α = Zoni-Shifted, β = 1/α  (Profile coefficients)\n";
        }
        else {
            std::cout << "Unknown profile coefficients\n";
        }

        std::cout << "------------------------------\n";

        std::cout << "nr_exp = " << nr_exp << ", nθ_exp = " << ntheta_exp << "\n";
        std::cout << "divideBy2 = " << divideBy2 << ", anisotropy = " << anisotropic_factor << "\n";
    }

    return true;
}

// Control Parameters
int ConfigParser::verbose() const
{
    return verbose_;
}
bool ConfigParser::paraview() const
{
    return paraview_;
}

int ConfigParser::maxOpenMPThreads() const
{
    return max_omp_threads_;
}
double ConfigParser::threadReductionFactor() const
{
    return thread_reduction_factor_;
}

bool ConfigParser::DirBC_Interior() const
{
    return DirBC_Interior_;
}

StencilDistributionMethod ConfigParser::stencilDistributionMethod() const
{
    return stencil_distribution_method_;
}
bool ConfigParser::cacheDensityProfileCoefficients() const
{
    return cache_density_profile_coefficients_;
}
bool ConfigParser::cacheDomainGeometry() const
{
    return cache_domain_geometry_;
}

const PolarGrid& ConfigParser::grid() const
{
    return grid_;
}

// Multigrid Parameters
bool ConfigParser::FMG() const
{
    return FMG_;
}
int ConfigParser::FMG_iterations() const
{
    return FMG_iterations_;
}
MultigridCycleType ConfigParser::FMG_cycle() const
{
    return FMG_cycle_;
}
ExtrapolationType ConfigParser::extrapolation() const
{
    return extrapolation_;
}
int ConfigParser::maxLevels() const
{
    return max_levels_;
}
int ConfigParser::preSmoothingSteps() const
{
    return pre_smoothing_steps_;
}
int ConfigParser::postSmoothingSteps() const
{
    return post_smoothing_steps_;
}
MultigridCycleType ConfigParser::multigridCycle() const
{
    return multigrid_cycle_;
}
int ConfigParser::maxIterations() const
{
    return max_iterations_;
}
ResidualNormType ConfigParser::residualNormType() const
{
    return residual_norm_type_;
}
std::optional<double> ConfigParser::absoluteTolerance() const
{
    return absolute_tolerance_;
}
std::optional<double> ConfigParser::relativeTolerance() const
{
    return relative_tolerance_;
}

const DomainGeometry& ConfigParser::domainGeometry() const
{
    return *domain_geometry_.get();
}

const DensityProfileCoefficients& ConfigParser::densityProfileCoefficients() const
{
    return *density_profile_coefficients_.get();
}

const BoundaryConditions& ConfigParser::boundaryConditions() const
{
    return *boundary_conditions_.get();
}

const SourceTerm& ConfigParser::sourceTerm() const
{
    return *source_term_.get();
}

const ExactSolution& ConfigParser::exactSolution() const
{
    return *exact_solution_.get();
}