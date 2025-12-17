
template <concepts::DomainGeometry DomainGeometry>
void GMGPolar<DomainGeometry>::setup()
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
        }
    }

    auto end_setup = std::chrono::high_resolution_clock::now();
    t_setup_total_ = std::chrono::duration<double>(end_setup - start_setup).count();
    LIKWID_STOP("Setup");
}
