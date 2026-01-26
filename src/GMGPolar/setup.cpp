#include "../../include/GMGPolar/gmgpolar.h"

int IGMGPolar::chooseNumberOfLevels(const PolarGrid& finestGrid)
{
    const int minRadialNodes      = 5;
    const int minAngularDivisions = 4;

    // Minimum level for Multigrid
    const int multigridMinLevel = 2;

    // Calculate radial maximum level
    int radialNodes    = finestGrid.nr();
    int radialMaxLevel = 1;
    while ((radialNodes + 1) / 2 >= minRadialNodes && (radialNodes + 1) % 2 == 0) {
        radialNodes = (radialNodes + 1) / 2;
        radialMaxLevel++;
    }

    // Calculate angular maximum level
    int angularDivisions = finestGrid.ntheta();
    int angularMaxLevel  = 1;
    while (angularDivisions / 2 >= minAngularDivisions && angularDivisions % 2 == 0 &&
           (angularDivisions / 2) % 2 == 0) {
        angularDivisions = angularDivisions / 2;
        angularMaxLevel++;
    }

    /* Currently unused: Number of levels which guarantee linear scalability */
    const int linear_complexity_levels = static_cast<int>(std::ceil(
        (2.0 * std::log(static_cast<double>(finestGrid.numberOfNodes())) - std::log(3.0)) / (3.0 * std::log(4.0))));

    // Determine the number of levels as the minimum of radial maximum level, angular maximum level,
    // and the maximum levels specified.
    int levels = std::min(radialMaxLevel, angularMaxLevel);
    if (max_levels_ > 0)
        levels = std::min(max_levels_, levels);

    // Check if levels is less than Multigrid minimum level and throw an error
    if (levels < multigridMinLevel) {
        throw std::runtime_error("Number of possible levels is less than Multigrid minimum level");
    }

    return levels;
}

void IGMGPolar::printSettings() const
{

    std::cout << "------------------------------\n";
    std::cout << "------- CMake Settings -------\n";
    std::cout << "------------------------------\n";
#ifdef NDEBUG
    std::cout << "Build: Release\n";
#else
    std::cout << "Build: Debug\n";
#endif

#ifdef GMGPOLAR_USE_MUMPS
    std::cout << "MUMPS: ON, ";
#else
    std::cout << "MUMPS: OFF, ";
#endif

#ifdef GMGPOLAR_USE_LIKWID
    std::cout << "Likwid: ON\n";
#else
    std::cout << "Likwid: OFF\n";
#endif

    std::cout << "------------------------------\n";
    std::cout << "------ General Settings ------\n";
    std::cout << "------------------------------\n";

    std::cout << "Maximum number of threads: " << max_omp_threads_ << "\n";

    if (DirBC_Interior_) {
        std::cout << "Dirichlet (Interior boundary condition)\n";
    }
    else {
        std::cout << "Across the origin (Interior boundary condition)\n";
    }

    if (stencil_distribution_method_ == StencilDistributionMethod::CPU_TAKE) {
        std::cout << "A-Take (Stencil Distribution)\n";
    }
    else {
        std::cout << "A-Give (Stencil Distribution)\n";
    }

    std::cout << "Domain geometry mode:" << " " << (cache_domain_geometry_ ? "Precomputed" : "On-the-fly") << "\n";

    std::cout << "Density profile mode:" << " " << (cache_density_profile_coefficients_ ? "Precomputed" : "On-the-fly")
              << "\n";

    std::cout << "------------------------------\n";
    std::cout << "---------- PolarGrid ---------\n";
    std::cout << "------------------------------\n";

    const PolarGrid& finest_grid   = levels_.front().grid();
    const PolarGrid& coarsest_grid = levels_.back().grid();

    std::cout << "r ∈ [" << finest_grid.radius(0) << ", " << finest_grid.radius(finest_grid.nr() - 1)
              << "], θ ∈ [0, 2π]\n";
    std::cout << "(nr × nθ) = (" << finest_grid.nr() << " × " << finest_grid.ntheta() << ") → (" << coarsest_grid.nr()
              << " × " << coarsest_grid.ntheta() << ")\n";
    std::cout << "Smoother: " << finest_grid.numberSmootherCircles() << " circles\n";
    std::cout << "Splitting radius = " << finest_grid.smootherSplittingRadius() << "\n";

    std::cout << "------------------------------\n";
    std::cout << "----- Multigrid settings -----\n";
    std::cout << "------------------------------\n";

    switch (extrapolation_) {
    case ExtrapolationType::NONE:
        std::cout << "No Extrapolation\n";
        break;
    case ExtrapolationType::IMPLICIT_EXTRAPOLATION:
        std::cout << "Implicit Extrapolation\n";
        break;
    case ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING:
        std::cout << "Implicit Extrapolation with Full Grid Smoothing\n";
        break;
    case ExtrapolationType::COMBINED:
        std::cout << "Combined Implicit Extrapolation\n";
        break;
    default:
        std::cout << "Unknown Extrapolation Type\n";
        break;
    }

    std::cout << "Multigrid Cycle: ";
    switch (multigrid_cycle_) {
    case MultigridCycleType::V_CYCLE:
        std::cout << "V(" << pre_smoothing_steps_ << "," << post_smoothing_steps_ << ")-Cycle\n";
        break;
    case MultigridCycleType::W_CYCLE:
        std::cout << "W(" << pre_smoothing_steps_ << "," << post_smoothing_steps_ << ")-Cycle\n";
        break;
    case MultigridCycleType::F_CYCLE:
        std::cout << "F(" << pre_smoothing_steps_ << "," << post_smoothing_steps_ << ")-Cycle\n";
        break;
    default:
        std::cout << "Unknown Cycle Type\n";
        break;
    }

    std::cout << "Number of levels: " << number_of_levels_ << "\n";

    std::cout << "Residual Norm: ";
    switch (residual_norm_type_) {
    case ResidualNormType::EUCLIDEAN:
        std::cout << "Euclidean (L2 norm)\n";
        break;
    case ResidualNormType::WEIGHTED_EUCLIDEAN:
        std::cout << "Weighted Euclidean (scaled L2 norm)\n";
        break;
    case ResidualNormType::INFINITY_NORM:
        std::cout << "Infinity Norm\n";
        break;
    default:
        std::cout << "Unknown Residual Norm Type\n";
        break;
    }

    std::cout << "Tolerances: abs = ";
    if (absolute_tolerance_) {
        std::cout << *absolute_tolerance_;
    }
    else {
        std::cout << "n/a";
    }
    std::cout << ", rel = ";
    if (relative_tolerance_) {
        std::cout << *relative_tolerance_;
    }
    else {
        std::cout << "n/a";
    }
    std::cout << ", maxiter = " << max_iterations_ << "\n";

    if (FMG_) {
        std::cout << "Full-Multigrid: " << FMG_iterations_ << "x ";
        switch (FMG_cycle_) {
        case MultigridCycleType::V_CYCLE:
            std::cout << "V(" << pre_smoothing_steps_ << "," << post_smoothing_steps_ << ")-Cycle\n";
            break;
        case MultigridCycleType::W_CYCLE:
            std::cout << "W(" << pre_smoothing_steps_ << "," << post_smoothing_steps_ << ")-Cycle\n";
            break;
        case MultigridCycleType::F_CYCLE:
            std::cout << "F(" << pre_smoothing_steps_ << "," << post_smoothing_steps_ << ")-Cycle\n";
            break;
        default:
            std::cout << "Unknown Configuration\n";
            break;
        }
    }
    else {
        std::cout << "Full-Multigrid: Disabled\n";
    }
}
