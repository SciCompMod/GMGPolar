#pragma once

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::setup()
{
    LIKWID_START("Setup");
    auto start_setup = std::chrono::high_resolution_clock::now();

    resetSetupPhaseTimings();

    auto start_setup_createLevels = std::chrono::high_resolution_clock::now();

    if (stencil_distribution_method_ == StencilDistributionMethod::TAKE) {
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
        writeToVTK("output_finest_grid", grid_);

    if (extrapolation_ != ExtrapolationType::NONE) {
        const double precision = 1e-12;
        if (!checkUniformRefinement(grid_, precision)) {
            std::cerr << "[Extrapolation Warning] Finest PolarGrid is not from a single uniform "
                         "refinement.\n";
        }
    }

    // ---------------------------------------------------------- //
    // Building PolarGrid and LevelCache for all multigrid levels //
    // ---------------------------------------------------------- //
    number_of_levels_ = chooseNumberOfLevels(grid_);
    levels_.clear();
    levels_.reserve(number_of_levels_);

    auto finest_levelCache = std::make_unique<LevelCache<DomainGeometry, DensityProfileCoefficients>>(
        *finest_grid, density_profile_coefficients_, domain_geometry_, cache_density_profile_coefficients_,
        cache_domain_geometry_);
    levels_.emplace_back(0, std::move(finest_grid), std::move(finest_levelCache), extrapolation_, FMG_, PCG_FMG_);

    for (int level_depth = 1; level_depth < number_of_levels_; level_depth++) {
        auto current_grid       = std::make_unique<PolarGrid>(coarseningGrid(levels_[level_depth - 1].grid()));
        auto current_levelCache = std::make_unique<LevelCache<DomainGeometry, DensityProfileCoefficients>>(
            *current_grid, density_profile_coefficients_, domain_geometry_, cache_density_profile_coefficients_,
            cache_domain_geometry_);
        levels_.emplace_back(level_depth, std::move(current_grid), std::move(current_levelCache), extrapolation_, FMG_,
                             PCG_FMG_);
    }

    auto end_setup_createLevels = std::chrono::high_resolution_clock::now();
    t_setup_createLevels_ = std::chrono::duration<double>(end_setup_createLevels - start_setup_createLevels).count();

    if (paraview_) {
        writeToVTK("output_coarsest_grid", levels_.back().grid());
    }

    // ------------------------------------- //
    // Initialize the interpolation operator //
    // ------------------------------------- //
    interpolation_ = std::make_unique<Interpolation>(DirBC_Interior_);

    if (verbose_ > 0) {
        printSettings(levels_[0].grid(), levels_[number_of_levels_ - 1].grid());
    }

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
        levels_[level_depth].initializeResidual(DirBC_Interior_, stencil_distribution_method_);
    }

    // ----------------------------------------- //
    // Build direct solver on the coarsest level //
    // ----------------------------------------- //
    auto start_setup_directSolver = std::chrono::high_resolution_clock::now();
    levels_[number_of_levels_ - 1].initializeDirectSolver(DirBC_Interior_, stencil_distribution_method_);
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
            levels_[0].initializeSmoothing(DirBC_Interior_, stencil_distribution_method_);
        }
        // PCG doesn't use extrapolated smoothing, so we only initialize it if PCG is disabled.
        if (do_extrapolated_smoothing && !PCG_) {
            levels_[0].initializeExtrapolatedSmoothing(DirBC_Interior_, stencil_distribution_method_);
        }
        for (int level_depth = 1; level_depth < number_of_levels_ - 1; level_depth++) {
            levels_[level_depth].initializeSmoothing(DirBC_Interior_, stencil_distribution_method_);
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
    constexpr int minRadialNodes      = 5;
    constexpr int minAngularDivisions = 4;

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

    // Determine the number of levels as the minimum of radial maximum level, angular maximum level,
    // and the maximum levels specified.
    int levels = std::min(radialMaxLevel, angularMaxLevel);
    if (max_levels_ > 0)
        levels = std::min(max_levels_, levels);

    // Extrapolation requires at least 2 levels
    if (extrapolation_ != ExtrapolationType::NONE && levels < 2) {
        std::cerr << "[GMGPolar Warning] Extrapolation disabled: requires at least 2 multigrid levels, but only "
                  << levels << " available.\n";

        extrapolation_ = ExtrapolationType::NONE;
    }

    return levels;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::discretize_rhs_f(
    const Level<DomainGeometry, DensityProfileCoefficients>& level, Vector<double> rhs_f)
{
    const PolarGrid& grid = level.grid();
    assert(std::ssize(rhs_f) == grid.numberOfNodes());

    const bool DirBC_Interior = DirBC_Interior_;

    if (level.levelCache().cacheDomainGeometry()) {
        /* DomainGeometry is cached */
        const ConstVector<double>& detDF_cache = level.levelCache().detDF();

        // ---------------------------------------------- //
        // Discretize rhs values (circular index section) //
        // ---------------------------------------------- //
        Kokkos::parallel_for(
            "discretize_rhs_f: Circular (cached)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>(
                {0, 0}, {grid.numberSmootherCircles(), grid.ntheta()}),
            KOKKOS_LAMBDA(const int i_r, const int i_theta) {
                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior)) {
                    const double h1    = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                    const double h2    = grid.radialSpacing(i_r);
                    const double k1    = grid.angularSpacing(i_theta - 1);
                    const double k2    = grid.angularSpacing(i_theta);
                    const double detDF = detDF_cache[grid.index(i_r, i_theta)];
                    rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * std::fabs(detDF);
                }
                else if (i_r == 0 && DirBC_Interior) {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
                else if (i_r == grid.nr() - 1) {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
            });

        // -------------------------------------------- //
        // Discretize rhs values (radial index section) //
        // -------------------------------------------- //
        Kokkos::parallel_for(
            "discretize_rhs_f: Radial (cached)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, grid.numberSmootherCircles()},
                                                                                  {grid.ntheta(), grid.nr()}),
            KOKKOS_LAMBDA(const int i_theta, const int i_r) {
                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior)) {
                    const double h1    = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                    const double h2    = grid.radialSpacing(i_r);
                    const double k1    = grid.angularSpacing(i_theta - 1);
                    const double k2    = grid.angularSpacing(i_theta);
                    const double detDF = detDF_cache[grid.index(i_r, i_theta)];
                    rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * std::fabs(detDF);
                }
                else if (i_r == 0 && DirBC_Interior) {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
                else if (i_r == grid.nr() - 1) {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
            });
    }
    else {
        /* DomainGeometry is not cached */
        // Local copy is required to avoid copying the class
        const DomainGeometry& domain_geometry = domain_geometry_;

        // ---------------------------------------------- //
        // Discretize rhs values (circular index section) //
        // ---------------------------------------------- //
        Kokkos::parallel_for(
            "discretize_rhs_f: Circular (uncached)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>(
                {0, 0}, {grid.numberSmootherCircles(), grid.ntheta()}),
            KOKKOS_LAMBDA(const int i_r, const int i_theta) {
                const double radius = grid.radius(i_r);
                const double theta  = grid.theta(i_theta);
                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior)) {
                    const double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                    const double h2 = grid.radialSpacing(i_r);
                    const double k1 = grid.angularSpacing(i_theta - 1);
                    const double k2 = grid.angularSpacing(i_theta);
                    /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                    /* The Jacobian matrix is: */
                    /* [Jrr, Jrt] */
                    /* [Jtr, Jtt] */
                    const double Jrr = domain_geometry.dFx_dr(radius, theta);
                    const double Jtr = domain_geometry.dFy_dr(radius, theta);
                    const double Jrt = domain_geometry.dFx_dtheta(radius, theta);
                    const double Jtt = domain_geometry.dFy_dtheta(radius, theta);
                    /* Compute the determinant of the Jacobian matrix */
                    const double detDF = Jrr * Jtt - Jrt * Jtr;
                    rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * std::fabs(detDF);
                }
                else if (i_r == 0 && DirBC_Interior) {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
                else if (i_r == grid.nr() - 1) {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
            });

        // -------------------------------------------- //
        // Discretize rhs values (radial index section) //
        // -------------------------------------------- //
        Kokkos::parallel_for(
            "discretize_rhs_f: Radial (uncached)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, grid.numberSmootherCircles()},
                                                                                  {grid.ntheta(), grid.nr()}),
            KOKKOS_LAMBDA(const int i_theta, const int i_r) {
                const double radius = grid.radius(i_r);
                const double theta  = grid.theta(i_theta);
                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior)) {
                    const double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                    const double h2 = grid.radialSpacing(i_r);
                    const double k1 = grid.angularSpacing(i_theta - 1);
                    const double k2 = grid.angularSpacing(i_theta);
                    /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                    /* The Jacobian matrix is: */
                    /* [Jrr, Jrt] */
                    /* [Jtr, Jtt] */
                    const double Jrr = domain_geometry.dFx_dr(radius, theta);
                    const double Jtr = domain_geometry.dFy_dr(radius, theta);
                    const double Jrt = domain_geometry.dFx_dtheta(radius, theta);
                    const double Jtt = domain_geometry.dFy_dtheta(radius, theta);
                    /* Compute the determinant of the Jacobian matrix */
                    const double detDF = Jrr * Jtt - Jrt * Jtr;
                    rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * std::fabs(detDF);
                }
                else if (i_r == 0 && DirBC_Interior) {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
                else if (i_r == grid.nr() - 1) {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
            });
    }

    Kokkos::fence();
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
template <concepts::BoundaryConditions BoundaryConditions, concepts::SourceTerm SourceTerm>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::build_rhs_f(
    const Level<DomainGeometry, DensityProfileCoefficients>& level, Vector<double> rhs_f,
    const BoundaryConditions& boundary_conditions, const SourceTerm& source_term)
{
    const PolarGrid& grid(level.grid());
    assert(std::ssize(rhs_f) == grid.numberOfNodes());

    const bool DirBC_Interior = DirBC_Interior_;

    // ----------------------------------------- //
    // Store rhs values (circular index section) //
    // ----------------------------------------- //
    Kokkos::parallel_for(
        "build_rhs_f: Circular",
        Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>(
            {0, 0}, {grid.numberSmootherCircles(), grid.ntheta()}),
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            const double radius = grid.radius(i_r);
            const double theta  = grid.theta(i_theta);
            if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior)) {
                rhs_f[grid.index(i_r, i_theta)] = source_term(i_r, i_theta);
            }
            else if (i_r == 0 && DirBC_Interior) {
                rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D_Interior(radius, theta);
            }
            else if (i_r == grid.nr() - 1) {
                rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D(radius, theta);
            }
        });

    // --------------------------------------- //
    // Store rhs values (radial index section) //
    // --------------------------------------- //
    Kokkos::parallel_for(
        "build_rhs_f: Radial",
        Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, grid.numberSmootherCircles()},
                                                                              {grid.ntheta(), grid.nr()}),
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            const double radius = grid.radius(i_r);
            const double theta  = grid.theta(i_theta);
            if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior)) {
                rhs_f[grid.index(i_r, i_theta)] = source_term(i_r, i_theta);
            }
            else if (i_r == 0 && DirBC_Interior) {
                rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D_Interior(radius, theta);
            }
            else if (i_r == grid.nr() - 1) {
                rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D(radius, theta);
            }
        });

    Kokkos::fence();
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

    std::cout << "Maximum number of threads: " << Kokkos::num_threads() << "\n";

    if (DirBC_Interior_) {
        std::cout << "Dirichlet (Interior boundary condition)\n";
    }
    else {
        std::cout << "Across the origin (Interior boundary condition)\n";
    }

    if (stencil_distribution_method_ == StencilDistributionMethod::TAKE) {
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

    HostConstVector<double> h_radius_finest = finest_grid.host_radii();
    HostConstVector<double> h_theta_finest  = finest_grid.host_theta();

    std::cout << "r ∈ [" << h_radius_finest(0) << ", " << h_radius_finest(finest_grid.nr() - 1) << "], θ ∈ [0, 2π]\n";
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

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
bool GMGPolar<DomainGeometry, DensityProfileCoefficients>::checkUniformRefinement(const PolarGrid& grid,
                                                                                  double tolerance) const
{
    HostConstVector<double> h_radius = grid.host_radii();
    HostConstVector<double> h_theta  = grid.host_theta();
    // Radial direction
    for (int i_r = 1; i_r < grid.nr() - 1; i_r += 2) {
        double left         = h_radius(i_r - 1);
        double right        = h_radius(i_r + 1);
        double expected_mid = 0.5 * (left + right);
        double actual_mid   = h_radius(i_r);

        double diff = std::abs(expected_mid - actual_mid);
        if (diff > tolerance) {
            std::cerr << "[Extrapolation Warning] Radial mismatch at i_r = " << i_r << "\n"
                      << "  left = " << left << ", right = " << right << "\n"
                      << "  expected = " << expected_mid << ", actual = " << actual_mid << "\n"
                      << "  diff = " << diff << " (tol = " << tolerance << ")\n";
            return false;
        }
    }

    // Angular direction
    for (int i_theta = 1; i_theta < grid.ntheta(); i_theta += 2) {
        double left         = h_theta(i_theta - 1);
        double right        = h_theta(i_theta + 1);
        double expected_mid = 0.5 * (left + right);
        double actual_mid   = h_theta(i_theta);

        double diff = std::abs(expected_mid - actual_mid);
        if (diff > tolerance) {
            std::cerr << "[Extrapolation Warning] Angular mismatch at i_theta = " << i_theta << "\n"
                      << "  left = " << left << ", right = " << right << "\n"
                      << "  expected = " << expected_mid << ", actual = " << actual_mid << "\n"
                      << "  diff = " << diff << " (tol = " << tolerance << ")\n";
            return false;
        }
    }

    return true;
}
