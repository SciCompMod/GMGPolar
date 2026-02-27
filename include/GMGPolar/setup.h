template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::setup()
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

    auto finest_levelCache = std::make_unique<LevelCache>(*finest_grid, density_profile_coefficients_, domain_geometry_,
                                                          cache_density_profile_coefficients_, cache_domain_geometry_);
    levels_.emplace_back(0, std::move(finest_grid), std::move(finest_levelCache), extrapolation_, FMG_);

    for (int level_depth = 1; level_depth < number_of_levels_; level_depth++) {
        auto current_grid = std::make_unique<PolarGrid>(coarseningGrid(levels_[level_depth - 1].grid()));
        auto current_levelCache =
            std::make_unique<LevelCache>(*current_grid, density_profile_coefficients_, domain_geometry_,
                                         cache_density_profile_coefficients_, cache_domain_geometry_);
        levels_.emplace_back(level_depth, std::move(current_grid), std::move(current_levelCache), extrapolation_, FMG_);
    }

    auto end_setup_createLevels = std::chrono::high_resolution_clock::now();
    t_setup_createLevels_ = std::chrono::duration<double>(end_setup_createLevels - start_setup_createLevels).count();

    if (paraview_)
        writeToVTK("output_coarsest_grid", levels_.back().grid());

    // ------------------------------------- //
    // Initialize the interpolation operator //
    // ------------------------------------- //
    interpolation_ = std::make_unique<Interpolation>(max_omp_threads_, DirBC_Interior_);

    if (verbose_ > 0)
        printSettings();

    // ------------------------------------------------ //
    // Define residual operator on all multigrid levels //
    // ------------------------------------------------ //
    for (int level_depth = 0; level_depth < number_of_levels_; level_depth++) {
        levels_[level_depth].initializeResidual(domain_geometry_, density_profile_coefficients_, DirBC_Interior_,
                                                max_omp_threads_, stencil_distribution_method_);
    }

    // ----------------------------------------- //
    // Build direct solver on the coarsest level //
    // ----------------------------------------- //
    auto start_setup_directSolver = std::chrono::high_resolution_clock::now();
    levels_[number_of_levels_ - 1].initializeDirectSolver(domain_geometry_, density_profile_coefficients_,
                                                          DirBC_Interior_, max_omp_threads_,
                                                          stencil_distribution_method_);
    auto end_setup_directSolver = std::chrono::high_resolution_clock::now();
    t_setup_directSolver_ += std::chrono::duration<double>(end_setup_directSolver - start_setup_directSolver).count();

    // ---------------------------------------------------------- //
    // Build the full-grid smoother and the extrapolated smoother //
    // ---------------------------------------------------------- //
    auto start_setup_smoother = std::chrono::high_resolution_clock::now();

    bool do_full_grid_smoothing    = false;
    bool do_extrapolated_smoothing = false;

    switch (extrapolation_) {

    case ExtrapolationType::NONE:
        do_full_grid_smoothing = true;
        break;

    case ExtrapolationType::IMPLICIT_EXTRAPOLATION:
        do_extrapolated_smoothing = true;
        break;

    case ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING:
        do_full_grid_smoothing = true;
        break;

    case ExtrapolationType::COMBINED:
        do_full_grid_smoothing    = true;
        do_extrapolated_smoothing = true;
        break;
    }

    full_grid_smoothing_ = do_full_grid_smoothing;

    if (number_of_levels_ > 1) {
        if (do_full_grid_smoothing) {
            levels_[0].initializeSmoothing(domain_geometry_, density_profile_coefficients_, DirBC_Interior_,
                                           max_omp_threads_, stencil_distribution_method_);
        }
        if (do_extrapolated_smoothing) {
            levels_[0].initializeExtrapolatedSmoothing(domain_geometry_, density_profile_coefficients_, DirBC_Interior_,
                                                       max_omp_threads_, stencil_distribution_method_);
        }
        for (int level_depth = 1; level_depth < number_of_levels_ - 1; level_depth++) {
            levels_[level_depth].initializeSmoothing(domain_geometry_, density_profile_coefficients_, DirBC_Interior_,
                                                     max_omp_threads_, stencil_distribution_method_);
        }
    }
    auto end_setup_smoother = std::chrono::high_resolution_clock::now();
    t_setup_smoother_ += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();

    auto end_setup = std::chrono::high_resolution_clock::now();
    t_setup_total_ = std::chrono::duration<double>(end_setup - start_setup).count();
    LIKWID_STOP("Setup");
}
