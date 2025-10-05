#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::setup()
{
    LIKWID_START("Setup");
    auto start_setup = std::chrono::high_resolution_clock::now();

    resetSetupPhaseTimings();

    auto start_setup_createLevels = std::chrono::high_resolution_clock::now();

    if (stencil_distribution_method_ == StencilDistributionMethod::CPU_TAKE) {
        if (!cache_density_profile_coefficients_ || !cache_domain_geometry_) {
            throw std::runtime_error("Error: Caching must be enabled for both density profile coefficients and domain "
                                     "geometry in 'Take' implementation strategy.");
        }
    }

    // -------------------------------- //
    // Create the finest mesh (level 0) //
    // -------------------------------- //
    auto finest_grid = std::make_unique<PolarGrid>(grid_);
    if (paraview_)
        writeToVTK("output_finest_grid", *finest_grid);

    if (extrapolation_ != ExtrapolationType::NONE) {
        // Check if grid comes from a single uniform refinement
        bool is_uniform_refinement = true;

        for (int i_r = 1; i_r < finest_grid->nr() - 1; i_r += 2) {
            double mid = 0.5 * (finest_grid->radius(i_r - 1) + finest_grid->radius(i_r + 1));
            if (std::abs(mid - finest_grid->radius(i_r)) > 1e-12) {
                is_uniform_refinement = false;
                break;
            }
        }
        for (int i_theta = 1; i_theta < finest_grid->ntheta(); i_theta += 2) {
            double mid = 0.5 * (finest_grid->theta(i_theta - 1) + finest_grid->theta(i_theta + 1));
            if (std::abs(mid - finest_grid->theta(i_theta)) > 1e-12) {
                is_uniform_refinement = false;
                break;
            }
        }

        if (!is_uniform_refinement) {
            throw std::runtime_error(
                "Extrapolation Error: Finest PolarGrid does not originate from a single uniform refinement.");
        }
    }

    // ---------------------------------------------------------- //
    // Building PolarGrid and LevelCache for all multigrid levels //
    // ---------------------------------------------------------- //
    number_of_levels_ = chooseNumberOfLevels(*finest_grid); /* Implementation below */
    levels_.clear();
    levels_.reserve(number_of_levels_);

    int level_depth        = 0;
    auto finest_levelCache = std::make_unique<LevelCache>(*finest_grid, density_profile_coefficients_, domain_geometry_,
                                                          cache_density_profile_coefficients_, cache_domain_geometry_);
    levels_.emplace_back(level_depth, std::move(finest_grid), std::move(finest_levelCache), extrapolation_, FMG_);

    for (level_depth = 1; level_depth < number_of_levels_; level_depth++) {
        auto current_grid       = std::make_unique<PolarGrid>(coarseningGrid(levels_[level_depth - 1].grid()));
        auto current_levelCache = std::make_unique<LevelCache>(levels_[level_depth - 1], *current_grid);
        levels_.emplace_back(level_depth, std::move(current_grid), std::move(current_levelCache), extrapolation_, FMG_);
    }

    auto end_setup_createLevels = std::chrono::high_resolution_clock::now();
    t_setup_createLevels_ = std::chrono::duration<double>(end_setup_createLevels - start_setup_createLevels).count();

    if (paraview_)
        writeToVTK("output_coarsest_grid", levels_.back().grid());

    // ----------------------------------------------------------- //
    // Initializing the optimal number of threads for OpenMP tasks //
    // ----------------------------------------------------------- //
    threads_per_level_.resize(number_of_levels_, max_omp_threads_);
    for (int level_depth = 0; level_depth < number_of_levels_; level_depth++) {
        threads_per_level_[level_depth] = std::max(
            1,
            std::min(max_omp_threads_,
                     static_cast<int>(std::floor(max_omp_threads_ * std::pow(thread_reduction_factor_, level_depth)))));
    }

    interpolation_ = std::make_unique<Interpolation>(threads_per_level_, DirBC_Interior_);

    if (verbose_ > 0)
        printSettings();

    // -------------------------------------------------------
    // Initializing various operators based on the level index
    for (int level_depth = 0; level_depth < number_of_levels_; level_depth++) {
        // ---------------------- //
        // Level 0 (finest Level) //
        // ---------------------- //
        if (level_depth == 0) {
            auto start_setup_smoother = std::chrono::high_resolution_clock::now();
            switch (extrapolation_) {
            case ExtrapolationType::NONE:
                full_grid_smoothing_ = true;
                levels_[level_depth].initializeSmoothing(domain_geometry_, density_profile_coefficients_,
                                                         DirBC_Interior_, threads_per_level_[level_depth],
                                                         stencil_distribution_method_);
                break;
            case ExtrapolationType::IMPLICIT_EXTRAPOLATION:
                full_grid_smoothing_ = false;
                levels_[level_depth].initializeExtrapolatedSmoothing(domain_geometry_, density_profile_coefficients_,
                                                                     DirBC_Interior_, threads_per_level_[level_depth],
                                                                     stencil_distribution_method_);
                break;
            case ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING:
                full_grid_smoothing_ = true;
                levels_[level_depth].initializeSmoothing(domain_geometry_, density_profile_coefficients_,
                                                         DirBC_Interior_, threads_per_level_[level_depth],
                                                         stencil_distribution_method_);
                break;
            case ExtrapolationType::COMBINED:
                full_grid_smoothing_ = true;
                levels_[level_depth].initializeSmoothing(domain_geometry_, density_profile_coefficients_,
                                                         DirBC_Interior_, threads_per_level_[level_depth],
                                                         stencil_distribution_method_);
                levels_[level_depth].initializeExtrapolatedSmoothing(domain_geometry_, density_profile_coefficients_,
                                                                     DirBC_Interior_, threads_per_level_[level_depth],
                                                                     stencil_distribution_method_);
                break;
            default:
                full_grid_smoothing_ = false;
                levels_[level_depth].initializeSmoothing(domain_geometry_, density_profile_coefficients_,
                                                         DirBC_Interior_, threads_per_level_[level_depth],
                                                         stencil_distribution_method_);
                levels_[level_depth].initializeExtrapolatedSmoothing(domain_geometry_, density_profile_coefficients_,
                                                                     DirBC_Interior_, threads_per_level_[level_depth],
                                                                     stencil_distribution_method_);
                break;
            }
            auto end_setup_smoother = std::chrono::high_resolution_clock::now();
            t_setup_smoother_ += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();
            levels_[level_depth].initializeResidual(domain_geometry_, density_profile_coefficients_, DirBC_Interior_,
                                                    threads_per_level_[level_depth], stencil_distribution_method_);
            levels_[level_depth].initializeSystemOperator(domain_geometry_, density_profile_coefficients_,
                                                          DirBC_Interior_, threads_per_level_[level_depth],
                                                          stencil_distribution_method_);
        }
        // -------------------------- //
        // Level n-1 (coarsest Level) //
        // -------------------------- //
        else if (level_depth == number_of_levels_ - 1) {
            auto start_setup_directSolver = std::chrono::high_resolution_clock::now();
            levels_[level_depth].initializeDirectSolver(domain_geometry_, density_profile_coefficients_,
                                                        DirBC_Interior_, threads_per_level_[level_depth],
                                                        stencil_distribution_method_);
            auto end_setup_directSolver = std::chrono::high_resolution_clock::now();
            t_setup_directSolver_ +=
                std::chrono::duration<double>(end_setup_directSolver - start_setup_directSolver).count();
            levels_[level_depth].initializeResidual(domain_geometry_, density_profile_coefficients_, DirBC_Interior_,
                                                    threads_per_level_[level_depth], stencil_distribution_method_);
            levels_[level_depth].initializeSystemOperator(domain_geometry_, density_profile_coefficients_,
                                                          DirBC_Interior_, threads_per_level_[level_depth],
                                                          stencil_distribution_method_);
        }
        // ------------------- //
        // Intermediate levels //
        // ------------------- //
        else {
            auto start_setup_smoother = std::chrono::high_resolution_clock::now();
            levels_[level_depth].initializeSmoothing(domain_geometry_, density_profile_coefficients_, DirBC_Interior_,
                                                     threads_per_level_[level_depth], stencil_distribution_method_);
            auto end_setup_smoother = std::chrono::high_resolution_clock::now();
            t_setup_smoother_ += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();
            levels_[level_depth].initializeResidual(domain_geometry_, density_profile_coefficients_, DirBC_Interior_,
                                                    threads_per_level_[level_depth], stencil_distribution_method_);
            levels_[level_depth].initializeSystemOperator(domain_geometry_, density_profile_coefficients_,
                                                          DirBC_Interior_, threads_per_level_[level_depth],
                                                          stencil_distribution_method_);
        }
    }

    auto end_setup = std::chrono::high_resolution_clock::now();
    t_setup_total_ = std::chrono::duration<double>(end_setup - start_setup).count();
    LIKWID_STOP("Setup");
}

int GMGPolar::chooseNumberOfLevels(const PolarGrid& finestGrid)
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
    const int linear_complexity_levels = std::ceil(
        (2.0 * std::log(static_cast<double>(finestGrid.numberOfNodes())) - std::log(3.0)) / (3.0 * std::log(4.0)));

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

void GMGPolar::printSettings() const
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

    std::cout << "r ∈ [" << finest_grid.radii().front() << ", " << finest_grid.radii().back() << "], θ ∈ [0, 2π]\n";
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
