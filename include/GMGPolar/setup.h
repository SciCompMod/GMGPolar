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

    auto finest_levelCache =
        std::make_unique<LevelCache<DomainGeometry>>(*finest_grid, density_profile_coefficients_, domain_geometry_,
                                                     cache_density_profile_coefficients_, cache_domain_geometry_);
    levels_.emplace_back(0, std::move(finest_grid), std::move(finest_levelCache), extrapolation_, FMG_, PCG_FMG_);

    for (int level_depth = 1; level_depth < number_of_levels_; level_depth++) {
        auto current_grid = std::make_unique<PolarGrid>(coarseningGrid(levels_[level_depth - 1].grid()));
        auto current_levelCache =
            std::make_unique<LevelCache<DomainGeometry>>(*current_grid, density_profile_coefficients_, domain_geometry_,
                                                         cache_density_profile_coefficients_, cache_domain_geometry_);
        levels_.emplace_back(level_depth, std::move(current_grid), std::move(current_levelCache), extrapolation_, FMG_,
                             PCG_FMG_);
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
        printSettings(levels_[0].grid(), levels_[number_of_levels_ - 1].grid());

    // ------------------------------- //
    // PCG-specific vector allocations //
    // ------------------------------- //
    if (PCG_) {
        const int grid_size = levels_[0].grid().numberOfNodes();
        if (std::ssize(pcg_solution_) != grid_size) {
            pcg_solution_ = Vector<double>("pcg_solution", grid_size);
        }
        if (std::ssize(pcg_search_direction_) != grid_size) {
            pcg_search_direction_ = Vector<double>("pcg_search_direction", grid_size);
        }
    }

    // ------------------------------------------------ //
    // Define residual operator on all multigrid levels //
    // ------------------------------------------------ //
    for (int level_depth = 0; level_depth < number_of_levels_; level_depth++) {
        levels_[level_depth].initializeResidual(density_profile_coefficients_, DirBC_Interior_, max_omp_threads_,
                                                stencil_distribution_method_);
    }

    // ----------------------------------------- //
    // Build direct solver on the coarsest level //
    // ----------------------------------------- //
    auto start_setup_directSolver = std::chrono::high_resolution_clock::now();
    levels_[number_of_levels_ - 1].initializeDirectSolver(density_profile_coefficients_, DirBC_Interior_,
                                                          max_omp_threads_, stencil_distribution_method_);
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
        // PCG uses non-extrapolated smoothing on level 0, so we need to initialize it if PCG is enabled.
        if (do_full_grid_smoothing || (PCG_ && PCG_MG_iterations_ > 0)) {
            levels_[0].initializeSmoothing(density_profile_coefficients_, DirBC_Interior_, max_omp_threads_,
                                           stencil_distribution_method_);
        }
        // PCG doesn't use extrapolated smoothing, so we only initialize it if PCG is disabled.
        if (do_extrapolated_smoothing && !PCG_) {
            levels_[0].initializeExtrapolatedSmoothing(density_profile_coefficients_, DirBC_Interior_, max_omp_threads_,
                                                       stencil_distribution_method_);
        }
        for (int level_depth = 1; level_depth < number_of_levels_ - 1; level_depth++) {
            levels_[level_depth].initializeSmoothing(density_profile_coefficients_, DirBC_Interior_, max_omp_threads_,
                                                     stencil_distribution_method_);
        }
    }
    auto end_setup_smoother = std::chrono::high_resolution_clock::now();
    t_setup_smoother_ += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();

    auto end_setup = std::chrono::high_resolution_clock::now();
    t_setup_total_ = std::chrono::duration<double>(end_setup - start_setup).count();
    LIKWID_STOP("Setup");
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
int GMGPolar<DomainGeometry, DensityProfileCoefficients>::chooseNumberOfLevels(const PolarGrid& finestGrid)
{
    const int minRadialNodes      = 5;
    const int minAngularDivisions = 4;

    // Minimum level for Multigrid
    const int multigridMinLevel = (extrapolation_ == ExtrapolationType::NONE) ? 1 : 2;

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

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::discretize_rhs_f(const Level<DomainGeometry>& level,
                                                                            Vector<double> rhs_f)
{
    const PolarGrid& grid = level.grid();
    assert(std::ssize(rhs_f) == grid.numberOfNodes());

    if (level.levelCache().cacheDomainGeometry()) {
        /* DomainGeometry is cached */
        const auto& detDF_cache = level.levelCache().detDF();
#pragma omp parallel
        {
// ---------------------------------------------- //
// Discretize rhs values (circular index section) //
// ---------------------------------------------- //
#pragma omp for nowait
            for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
                double r = grid.radius(i_r);
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    double theta = grid.theta(i_theta);
                    if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                        double h1          = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                        double h2          = grid.radialSpacing(i_r);
                        double k1          = grid.angularSpacing(i_theta - 1);
                        double k2          = grid.angularSpacing(i_theta);
                        const double detDF = detDF_cache[grid.index(i_r, i_theta)];
                        rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
                    }
                    else if (i_r == 0 && DirBC_Interior_) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                    else if (i_r == grid.nr() - 1) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                }
            }

// -------------------------------------------- //
// Discretize rhs values (radial index section) //
// -------------------------------------------- //
#pragma omp for nowait
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                double theta = grid.theta(i_theta);
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    double r = grid.radius(i_r);
                    if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                        double h1          = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                        double h2          = grid.radialSpacing(i_r);
                        double k1          = grid.angularSpacing(i_theta - 1);
                        double k2          = grid.angularSpacing(i_theta);
                        const double detDF = detDF_cache[grid.index(i_r, i_theta)];
                        rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
                    }
                    else if (i_r == 0 && DirBC_Interior_) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                    else if (i_r == grid.nr() - 1) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                }
            }
        }
    }
    else {
        /* DomainGeometry is not cached */

#pragma omp parallel
        {
// ---------------------------------------------- //
// Discretize rhs values (circular index section) //
// ---------------------------------------------- //
#pragma omp for nowait
            for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
                double r = grid.radius(i_r);
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    double theta = grid.theta(i_theta);

                    if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                        double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                        double h2 = grid.radialSpacing(i_r);
                        double k1 = grid.angularSpacing(i_theta - 1);
                        double k2 = grid.angularSpacing(i_theta);
                        /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                        /* The Jacobian matrix is: */
                        /* [Jrr, Jrt] */
                        /* [Jtr, Jtt] */
                        double Jrr = domain_geometry_.dFx_dr(r, theta);
                        double Jtr = domain_geometry_.dFy_dr(r, theta);
                        double Jrt = domain_geometry_.dFx_dt(r, theta);
                        double Jtt = domain_geometry_.dFy_dt(r, theta);
                        /* Compute the determinant of the Jacobian matrix */
                        double detDF = Jrr * Jtt - Jrt * Jtr;
                        rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
                    }
                    else if (i_r == 0 && DirBC_Interior_) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                    else if (i_r == grid.nr() - 1) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                }
            }

// -------------------------------------------- //
// Discretize rhs values (radial index section) //
// -------------------------------------------- //
#pragma omp for nowait
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                double theta = grid.theta(i_theta);

                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    double r = grid.radius(i_r);
                    if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                        double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                        double h2 = grid.radialSpacing(i_r);
                        double k1 = grid.angularSpacing(i_theta - 1);
                        double k2 = grid.angularSpacing(i_theta);
                        /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                        /* The Jacobian matrix is: */
                        /* [Jrr, Jrt] */
                        /* [Jtr, Jtt] */
                        double Jrr = domain_geometry_.dFx_dr(r, theta);
                        double Jtr = domain_geometry_.dFy_dr(r, theta);
                        double Jrt = domain_geometry_.dFx_dt(r, theta);
                        double Jtt = domain_geometry_.dFy_dt(r, theta);
                        /* Compute the determinant of the Jacobian matrix */
                        double detDF = Jrr * Jtt - Jrt * Jtr;
                        rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
                    }
                    else if (i_r == 0 && DirBC_Interior_) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                    else if (i_r == grid.nr() - 1) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                }
            }
        }
    }
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
template <concepts::BoundaryConditions BoundaryConditions>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::build_rhs_f(const Level<DomainGeometry>& level,
                                                                       Vector<double> rhs_f,
                                                                       const BoundaryConditions& boundary_conditions,
                                                                       const SourceTerm& source_term)
{
    const PolarGrid& grid = level.grid();
    assert(std::ssize(rhs_f) == grid.numberOfNodes());

#pragma omp parallel
    {
// ----------------------------------------- //
// Store rhs values (circular index section) //
// ----------------------------------------- //
#pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                double theta = grid.theta(i_theta);

                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                    rhs_f[grid.index(i_r, i_theta)] = source_term(i_r, i_theta);
                }
                else if (i_r == 0 && DirBC_Interior_) {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D_Interior(r, theta);
                }
                else if (i_r == grid.nr() - 1) {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D(r, theta);
                }
            }
        }

// --------------------------------------- //
// Store rhs values (radial index section) //
// --------------------------------------- //
#pragma omp for
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            double theta = grid.theta(i_theta);

            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                double r = grid.radius(i_r);
                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                    rhs_f[grid.index(i_r, i_theta)] = source_term(i_r, i_theta);
                }
                else if (i_r == 0 && DirBC_Interior_) {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D_Interior(r, theta);
                }
                else if (i_r == grid.nr() - 1) {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D(r, theta);
                }
            }
        }
    }
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::printSettings(const PolarGrid& finest_grid,
                                                                         const PolarGrid& coarsest_grid) const
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

    if (PCG_) {
        std::cout << "Preconditioned Conjugate Gradient: Enabled\n";
        if (PCG_FMG_) {
            std::cout << "- PCG Full-Multigrid: " << PCG_FMG_iterations_ << "x ";
            switch (PCG_FMG_cycle_) {
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
            std::cout << "PCG Full-Multigrid: Disabled\n";
        }

        std::cout << "- PCG Multigrid Iteration: " << PCG_MG_iterations_ << "x ";
        switch (PCG_MG_cycle_) {
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
        std::cout << "Preconditioned Conjugate Gradient: Disabled\n";
    }
}
