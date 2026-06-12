#pragma once

// =============================================================================
//   Main Solver Routine
// =============================================================================
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
template <concepts::BoundaryConditions BoundaryConditions, concepts::SourceTerm SourceTerm>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::solve(const BoundaryConditions& boundary_conditions,
                                                                 const SourceTerm& source_term)
{
    auto start_setup_rhs = std::chrono::high_resolution_clock::now();

    /* ------------------------------------- */
    /* Build rhs_f on Level 0 (finest Level) */
    /* ------------------------------------- */
    Vector<double> rhs_f = levels_[0].rhs();
    build_rhs_f(levels_[0], rhs_f, boundary_conditions, source_term);

    /* ---------------- */
    /* Discretize rhs_f */
    /* ---------------- */
    int initial_rhs_f_levels = FMG_ ? number_of_levels_ : (extrapolation_ == ExtrapolationType::NONE ? 1 : 2);
    // Loop through the levels, injecting and discretizing rhs
    for (int level_depth = 0; level_depth < initial_rhs_f_levels; ++level_depth) {
        Level<DomainGeometry, DensityProfileCoefficients>& current_level = levels_[level_depth];
        // Inject rhs if there is a next level
        if (level_depth + 1 < initial_rhs_f_levels) {
            Level<DomainGeometry, DensityProfileCoefficients>& next_level = levels_[level_depth + 1];
            injection(level_depth, next_level.rhs(), current_level.rhs());
        }
        // Discretize the rhs for the current level
		auto h_rhs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, current_level.rhs());
        discretize_rhs_f(current_level, h_rhs);
		Kokkos::deep_copy(current_level.rhs(), h_rhs);
    }

    auto end_setup_rhs = std::chrono::high_resolution_clock::now();
    t_setup_rhs_       = std::chrono::duration<double>(end_setup_rhs - start_setup_rhs).count();

    LIKWID_START("Solve");
    auto start_solve = std::chrono::high_resolution_clock::now();

    // Clear solve-phase timings
    resetSolvePhaseTimings();

    /* ---------------------------- */
    /* Initialize starting solution */
    /* ---------------------------- */
    auto start_initial_approximation = std::chrono::high_resolution_clock::now();

    if (!FMG_) {
        // Assign zero initial guess if not using FMG
        assign(levels_[0].solution(), 0.0);
    }
    else {
        // Compute an initial approximation to the discrete system
        //
        //     A * levels_[0].solution() = levels_[0].rhs()
        //
        // using the Full Multigrid (FMG) algorithm.
        //
        // Prerequisite:
        // The right-hand side must already be properly initialized on all levels,
        // i.e., constructed on the finest level and transferred to coarser levels
        // via injection and discretization. This ensures consistency of the coarse-
        // grid operators and guarantees correctness of the multigrid hierarchy.
        //
        // The FMG algorithm performs an exact solve on the coarsest grid and then
        // successively prolongates the solution to finer grids, applying multigrid
        // cycles on each level to reduce the error. This produces a high-quality
        // initial guess on the finest level, typically reducing the error to the
        // order of the discretization error and significantly accelerating convergence
        // of the subsequent solve phase.
        fullMultigridApproximation(FMG_cycle_, FMG_iterations_);
    }

    auto end_initial_approximation = std::chrono::high_resolution_clock::now();
    t_solve_initial_approximation_ =
        std::chrono::duration<double>(end_initial_approximation - start_initial_approximation).count();

    // These times are included in the initial approximation and don't count towards the multigrid cyclces.
    resetAvgMultigridCycleTimings();

    /* --------------------------------------- */
    /* Start Solver at finest level (depth 0)  */
    /* --------------------------------------- */
    Level<DomainGeometry, DensityProfileCoefficients>& level = levels_[0];

    number_of_iterations_                 = 0;
    double initial_residual_norm          = 1.0;
    double current_residual_norm          = 1.0;
    double current_relative_residual_norm = 1.0;

    printIterationHeader(exact_solution_);

    /* ---------------------------------------------- */
    /* Test solution against exact solution if given. */
    /* ---------------------------------------------- */
    LIKWID_STOP("Solver");
    auto start_check_exact_error = std::chrono::high_resolution_clock::now();

    if (exact_solution_) {
        /* Fill analytical solution on the host. */
        // Note: We must compute the analytical solution values on the host because the
        //       ExactSolution object is not designed to be used on the device.
        //       As it is planned that the solution and error vector will use device memory during GPU porting,
        //       we need to transfer the exact solution values from host to device using the Kokkos::deep_copy operation.
        if (std::ssize(analytical_solution_host_) != level.grid().numberOfNodes()) {
            analytical_solution_host_ = HostVector<double>("Analytical Solution", level.solution().size());
        }
        const PolarGrid<Kokkos::HostSpace> grid(level.grid());
        computeAnalyticalSolutionOnHost(grid, analytical_solution_host_, *exact_solution_);

        // Evaluate the error of the initial approximation against the exact solution, if provided.
        // We use level.residual() as a temporary vector to store the error values.
        exact_errors_.push_back(
            evaluateExactError(level.grid(), level.solution(), analytical_solution_host_, level.residual()));
    }
    auto end_check_exact_error = std::chrono::high_resolution_clock::now();
    t_check_exact_error_ += std::chrono::duration<double>(end_check_exact_error - start_check_exact_error).count();
    LIKWID_START("Solver");

    /* ---------------------------- */
    /* Compute convergence criteria */
    /* ---------------------------- */
    auto start_check_convergence = std::chrono::high_resolution_clock::now();

    // Initializes level.residual() and sets up the convergence criteria.
    updateResidualNorms(level, number_of_iterations_, initial_residual_norm, current_residual_norm,
                        current_relative_residual_norm);

    auto end_check_convergence = std::chrono::high_resolution_clock::now();
    t_check_convergence_ += std::chrono::duration<double>(end_check_convergence - start_check_convergence).count();

    printIterationInfo(number_of_iterations_, current_residual_norm, current_relative_residual_norm, exact_solution_);

    if (!converged(current_residual_norm, current_relative_residual_norm)) {
        if (!PCG_) {
            // Solve A*x = b directly using multigrid cycles (V/W/F-cycle)
            // until convergence or max_iterations_ is reached.
            solveMultigrid(initial_residual_norm, current_residual_norm, current_relative_residual_norm);
        }
        else {
            // Solve A*x = b using Preconditioned Conjugate Gradient (PCG),
            // with multigrid cycles as the preconditioner (i.e., one multigrid
            // cycle approximates the action of A^{-1} at each PCG iteration).
            auto start_conjugate_gradient = std::chrono::high_resolution_clock::now();

            solvePCG(initial_residual_norm, current_residual_norm, current_relative_residual_norm);

            auto end_conjugate_gradient = std::chrono::high_resolution_clock::now();
            t_conjugate_gradient_ +=
                std::chrono::duration<double>(end_conjugate_gradient - start_conjugate_gradient).count();
        }
    }

    /* ---------------------- */
    /* Post-solution analysis */
    /* ---------------------- */
    if (number_of_iterations_ > 0) {
        /* --------------------------------------------- */
        /* Compute the average Multigrid Iteration times */
        /* --------------------------------------------- */
        t_avg_MGC_total_ = t_solve_multigrid_iterations_ / number_of_iterations_;
        t_avg_MGC_preSmoothing_ /= number_of_iterations_;
        t_avg_MGC_postSmoothing_ /= number_of_iterations_;
        t_avg_MGC_residual_ /= number_of_iterations_;
        t_avg_MGC_directSolver_ /= number_of_iterations_;

        /* -------------------------------- */
        /* Compute the reduction factor rho */
        /* -------------------------------- */
        mean_residual_reduction_factor_ =
            std::pow(current_residual_norm / initial_residual_norm, 1.0 / number_of_iterations_);

        if (verbose_ > 0) {
            std::cout << "------------------------------\n";
            std::cout << "Total Iterations: " << number_of_iterations_ << "\n";
            std::cout << "Reduction Factor: ρ = " << mean_residual_reduction_factor_ << "\n";
        }
    }

    auto end_solve = std::chrono::high_resolution_clock::now();
    t_solve_total_ = std::chrono::duration<double>(end_solve - start_solve).count() - t_check_exact_error_;
    LIKWID_STOP("Solve");

    if (paraview_) {
        writeToVTK("output_solution", level, level.solution());
        if (exact_solution_) {
            evaluateExactError(level.grid(), level.solution(), analytical_solution_host_, level.residual());
            writeToVTK("output_error", level, level.residual());
        }
    }
}

// =============================================================================
//   Full Multigrid Approximation
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::fullMultigridApproximation(MultigridCycleType FMG_cycle,
                                                                                      int FMG_iterations)
{
    // Start from the coarsest level
    int coarsest_depth                                                = number_of_levels_ - 1;
    Level<DomainGeometry, DensityProfileCoefficients>& coarsest_level = levels_[coarsest_depth];

    // Solve directly on the coarsest level
    Kokkos::deep_copy(coarsest_level.solution(), coarsest_level.rhs());
	auto coarsest_level_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), coarsest_level.solution());
    coarsest_level.directSolveInPlace(coarsest_level_solution);
	Kokkos::deep_copy(coarsest_level.solution(), coarsest_level_solution);

    // Prolongate the solution from the coarsest level up to the finest, while applying Multigrid Cycles on each level
    for (int depth = coarsest_depth; depth > 0; --depth) {
        Level<DomainGeometry, DensityProfileCoefficients>& coarse_level = levels_[depth]; // Current coarse level
        Level<DomainGeometry, DensityProfileCoefficients>& fine_level   = levels_[depth - 1]; // Next finer level

        // The bi-cubic FMG interpolation is of higher order
        FMGInterpolation(coarse_level.level_depth(), fine_level.solution(), coarse_level.solution());

        // Apply some FMG iterations, except on the finest level,
        // where the interpolated solution is sufficiently accurate as an initial guess
        if (fine_level.level_depth() > 0) {
            applyMultigridIterations(fine_level, FMG_cycle, FMG_iterations);
        }
    }
}

// =============================================================================
//   Solver Loops
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::solveMultigrid(double& initial_residual_norm,
                                                                          double& current_residual_norm,
                                                                          double& current_relative_residual_norm)
{
    Level<DomainGeometry, DensityProfileCoefficients>& level = levels_[0];

    while (number_of_iterations_ < max_iterations_) {
        /* ----------------------- */
        /* Perform Multigrid Cycle */
        /* ----------------------- */
        auto start_solve_multigrid_iterations = std::chrono::high_resolution_clock::now();

        if (extrapolation_ == ExtrapolationType::NONE) {
            applyMultigridIterations(level, multigrid_cycle_, 1);
        }
        else {
            applyExtrapolatedMultigridIterations(level, multigrid_cycle_, 1);
        }

        number_of_iterations_++;

        auto end_solve_multigrid_iterations = std::chrono::high_resolution_clock::now();
        t_solve_multigrid_iterations_ +=
            std::chrono::duration<double>(end_solve_multigrid_iterations - start_solve_multigrid_iterations).count();

        /* ---------------------------------------------- */
        /* Test solution against exact solution if given. */
        /* ---------------------------------------------- */
        LIKWID_STOP("Solver");
        auto start_check_exact_error = std::chrono::high_resolution_clock::now();

        if (exact_solution_) {
            exact_errors_.push_back(
                evaluateExactError(level.grid(), level.solution(), analytical_solution_host_, level.residual()));
        }

        auto end_check_exact_error = std::chrono::high_resolution_clock::now();
        t_check_exact_error_ += std::chrono::duration<double>(end_check_exact_error - start_check_exact_error).count();
        LIKWID_START("Solver");

        /* ---------------------------- */
        /* Compute convergence criteria */
        /* ---------------------------- */
        auto start_check_convergence = std::chrono::high_resolution_clock::now();

        if (absolute_tolerance_.has_value() || relative_tolerance_.has_value()) {
            updateResidualNorms(level, number_of_iterations_, initial_residual_norm, current_residual_norm,
                                current_relative_residual_norm);
        }

        auto end_check_convergence = std::chrono::high_resolution_clock::now();
        t_check_convergence_ += std::chrono::duration<double>(end_check_convergence - start_check_convergence).count();

        printIterationInfo(number_of_iterations_, current_residual_norm, current_relative_residual_norm,
                           exact_solution_);

        if (converged(current_residual_norm, current_relative_residual_norm))
            break;
    }
}

/* -------------------------------------------------- */
/* Vectors for PCG (Preconditioned Conjugate Gradient)
 * https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
 *
 * Dedicated vectors:
 *   x  (solution)            -> pcg_solution_
 *   p  (search direction)    -> pcg_search_direction_
 *
 * Reused vectors (to avoid extra allocations):
 *   r    (residual)                       -> level.rhs()
 *   z    (preconditioned residual)        -> level.solution()
 *   A*p  (matrix applied to search dir.)  -> level.residual()
 */
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::solvePCG(double& initial_residual_norm,
                                                                    double& current_residual_norm,
                                                                    double& current_relative_residual_norm)
{
    Level<DomainGeometry, DensityProfileCoefficients>& level = levels_[0];

    // x = initial guess
    Kokkos::deep_copy(pcg_solution_, level.solution());
    // r = residual of initial guess
    Kokkos::deep_copy(level.rhs(), level.residual());

    // z = M^{-1} * r (preconditioned residual)
    if (PCG_FMG_) {
		auto h_rhs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, level.rhs());
        initRhsHierarchy(h_rhs);
		Kokkos::deep_copy(level.rhs(), h_rhs);
        fullMultigridApproximation(PCG_FMG_cycle_, PCG_FMG_iterations_);
    }
    else {
        // z = I^{-1} * r (no preconditioning)
        Kokkos::deep_copy(level.solution(), level.rhs());
    }
    applyMultigridIterations(level, PCG_MG_cycle_, PCG_MG_iterations_);

    // p = z
    Kokkos::deep_copy(pcg_search_direction_, level.solution());

    // r^T * z
	auto h_rhs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, level.rhs());
    double r_z = dot_product(HostConstVector<double>(h_rhs), HostConstVector<double>(level.solution()));

    while (number_of_iterations_ < max_iterations_) {

        // A_p = A * p
        auto level_residual       = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), level.residual());
        auto pcg_search_direction = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), pcg_search_direction_);
        level.applySystemOperator(level_residual, pcg_search_direction);
        Kokkos::deep_copy(level.residual(), level_residual);
        if (extrapolation_ != ExtrapolationType::NONE) {
            assert(number_of_levels_ > 1);
            Level<DomainGeometry, DensityProfileCoefficients>& next_level = levels_[level.level_depth() + 1];
            auto next_level_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), next_level.solution());
            auto next_level_residual = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), next_level.residual());
            injection(0, next_level_solution, pcg_search_direction);
            next_level.applySystemOperator(next_level_residual, next_level_solution);
            applyExtrapolation(0, level_residual, next_level_residual);
            Kokkos::deep_copy(level.residual(), level_residual);
            Kokkos::deep_copy(next_level.residual(), next_level_residual);
        }

        // alpha = (r^T * z) / (p^T * A*p)
        double alpha = r_z / dot_product(HostConstVector<double>(pcg_search_direction_),
                                         HostConstVector<double>(level.residual()));

        // x += alpha * p
        linear_combination(pcg_solution_, 1.0, HostConstVector<double>(pcg_search_direction_), alpha);

        // r -= alpha * A*p
		Kokkos::deep_copy(h_rhs, level.rhs());
        linear_combination(h_rhs, 1.0, HostConstVector<double>(level.residual()), -alpha);

        /* ---------------------------- */
        /* Compute convergence criteria */
        /* ---------------------------- */
        auto start_check_convergence = std::chrono::high_resolution_clock::now();

        current_residual_norm = residualNorm(residual_norm_type_, level, h_rhs);
        residual_norms_.push_back(current_residual_norm);
        current_relative_residual_norm = current_residual_norm / initial_residual_norm;

        auto end_check_convergence = std::chrono::high_resolution_clock::now();
        t_check_convergence_ += std::chrono::duration<double>(end_check_convergence - start_check_convergence).count();

        /* ---------------------------------------------- */
        /* Test solution against exact solution if given. */
        /* ---------------------------------------------- */
        LIKWID_STOP("Solver");
        auto start_check_exact_error = std::chrono::high_resolution_clock::now();

        if (exact_solution_) {
            exact_errors_.push_back(
                evaluateExactError(level.grid(), pcg_solution_, analytical_solution_host_, level.residual()));
        }

        auto end_check_exact_error = std::chrono::high_resolution_clock::now();
        t_check_exact_error_ += std::chrono::duration<double>(end_check_exact_error - start_check_exact_error).count();
        LIKWID_START("Solver");

        number_of_iterations_++;

        printIterationInfo(number_of_iterations_, current_residual_norm, current_relative_residual_norm,
                           exact_solution_);

        if (converged(current_residual_norm, current_relative_residual_norm))
            break;

        // z = M^{-1} * r (preconditioned residual)
        if (PCG_FMG_) {
            initRhsHierarchy(h_rhs);
            fullMultigridApproximation(PCG_FMG_cycle_, PCG_FMG_iterations_);
        }
        else {
            // z = I^{-1} * r (no preconditioning)
            Kokkos::deep_copy(level.solution(), h_rhs);
        }
        applyMultigridIterations(level, PCG_MG_cycle_, PCG_MG_iterations_);

        // r^T * z
        double r_z_new = dot_product(HostConstVector<double>(h_rhs), HostConstVector<double>(level.solution()));
        // beta = (current r^T * z) / (previous r^T * z)
        double beta = r_z_new / r_z;
        // r_z = r^T * z for next iteration
        r_z = r_z_new;

        // p *= beta
        multiply(pcg_search_direction_, beta);
        // p += z
        add(pcg_search_direction_, HostConstVector<double>(level.solution()));
    }
    // level.solution() = x
    Kokkos::deep_copy(level.solution(), pcg_solution_);
}

// =============================================================================
//   Apply Multigrid Iterations
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::applyMultigridIterations(
    Level<DomainGeometry, DensityProfileCoefficients>& level, MultigridCycleType cycle, int iterations)
{
    for (int i = 0; i < iterations; i++) {
	    auto h_rhs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, level.rhs());
		auto residual = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace{}, level.residual());
		auto solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace{}, level.solution());
        switch (cycle) {
        case MultigridCycleType::V_CYCLE:
            multigrid_V_Cycle(level.level_depth(), level.solution(), h_rhs, level.residual());
            break;
        case MultigridCycleType::W_CYCLE:
            multigrid_W_Cycle(level.level_depth(), solution, level.rhs(), residual);
			Kokkos::deep_copy(level.residual(), residual);
			Kokkos::deep_copy(level.solution(), solution);
            break;
        case MultigridCycleType::F_CYCLE:
            multigrid_F_Cycle(level.level_depth(), level.solution(), h_rhs, level.residual());
            break;
        default:
            std::cerr << "Error: Unknown multigrid cycle type!" << std::endl;
            break;
        }
    }
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::applyExtrapolatedMultigridIterations(
    Level<DomainGeometry, DensityProfileCoefficients>& level, MultigridCycleType cycle, int iterations)
{
    for (int i = 0; i < iterations; i++) {
		auto h_rhs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, level.rhs());
		auto residual = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace{}, level.residual());
		auto solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace{}, level.solution());
        switch (cycle) {
        case MultigridCycleType::V_CYCLE:
            extrapolated_multigrid_V_Cycle(level.level_depth(), solution, level.rhs(), residual);
			Kokkos::deep_copy(level.residual(), residual);
			Kokkos::deep_copy(level.solution(), solution);
            break;
        case MultigridCycleType::W_CYCLE:
            extrapolated_multigrid_W_Cycle(level.level_depth(), solution, level.rhs(), residual);
			Kokkos::deep_copy(level.residual(), residual);
			Kokkos::deep_copy(level.solution(), solution);
            break;
        case MultigridCycleType::F_CYCLE:
            extrapolated_multigrid_F_Cycle(level.level_depth(), level.solution(), h_rhs, level.residual());
            break;
        default:
            std::cerr << "Error: Unknown multigrid cycle type!" << std::endl;
            break;
        }
    }
}

// =============================================================================
//   Residual Handling Functions
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::updateResidualNorms(
    Level<DomainGeometry, DensityProfileCoefficients>& level, int iteration, double& initial_residual_norm,
    double& current_residual_norm, double& current_relative_residual_norm)
{
    auto solution            = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), level.solution());
    auto residual      = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), level.residual());
    level.computeResidual(residual, level.rhs(), solution);
    if (extrapolation_ != ExtrapolationType::NONE) {
        Level<DomainGeometry, DensityProfileCoefficients>& next_level = levels_[level.level_depth() + 1];
        auto next_level_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), next_level.solution());
        auto next_level_residual = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), next_level.residual());
        auto next_level_rhs = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), next_level.rhs());
        injection(level.level_depth(), next_level_solution, solution);
        next_level.computeResidual(next_level_residual, next_level_rhs, next_level_solution);
        applyExtrapolation(level.level_depth(), residual, next_level_residual);
        Kokkos::deep_copy(next_level.solution(), next_level_solution);
    }
    Kokkos::deep_copy(level.residual(), residual);
    Kokkos::deep_copy(level.solution(), solution);

    current_residual_norm = residualNorm(residual_norm_type_, level, level.residual());
    residual_norms_.push_back(current_residual_norm);

    if (number_of_iterations_ == 0) {
	    auto h_rhs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, level.rhs());
        initial_residual_norm = !FMG_ ? current_residual_norm : residualNorm(residual_norm_type_, level, h_rhs);
    }
    current_relative_residual_norm = current_residual_norm / initial_residual_norm;

    // Combined Smoothing: If small residual reduction, turn off full grid smoothing.
    if (number_of_iterations_ > 0) {
        const double convergence_factor = 0.7;
        const double current_residual_reduction_factor =
            residual_norms_[number_of_iterations_] / residual_norms_[number_of_iterations_ - 1];
        if (current_residual_reduction_factor > convergence_factor && extrapolation_ == ExtrapolationType::COMBINED &&
            full_grid_smoothing_) {
            full_grid_smoothing_ = false;
            if (verbose_ > 0)
                std::cout << "Switching from full grid smoothing to standard extrapolated smoothing.\n";
        }
    }
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
double GMGPolar<DomainGeometry, DensityProfileCoefficients>::residualNorm(
    const ResidualNormType& norm_type, const Level<DomainGeometry, DensityProfileCoefficients>& level,
    HostConstVector<double> residual) const
{
    switch (norm_type) {
    case ResidualNormType::EUCLIDEAN:
        return l2_norm(residual);
    case ResidualNormType::WEIGHTED_EUCLIDEAN:
        return l2_norm(residual) / std::sqrt(level.grid().numberOfNodes());
    case ResidualNormType::INFINITY_NORM:
        return infinity_norm(residual);
    default:
        throw std::invalid_argument("Unknown ResidualNormType");
    }
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::applyExtrapolation(int current_level,
                                                                              Vector<double> fine_values,
                                                                              ConstVector<double> coarse_values)
{
    const PolarGrid<DefaultMemorySpace>& fineGrid(levels_[current_level].grid());
    const PolarGrid<DefaultMemorySpace>& coarseGrid(levels_[current_level + 1].grid());

    assert(std::ssize(fine_values) == fineGrid.numberOfNodes());
    assert(std::ssize(coarse_values) == coarseGrid.numberOfNodes());

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    /* Circular Indexing Section */
    /* For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Extrapolation: Apply Extrapolation (Circular)",
        Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>(
            {0, 0}, {fineGrid.numberSmootherCircles(), fineGrid.ntheta()}),
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            const int fine_idx = fineGrid.index(i_r, i_theta);
            if (i_r & 1 || i_theta & 1) {
                fine_values[fine_idx] *= 4.0 / 3.0;
            }
            else {
                const int coarse_idx  = coarseGrid.index(i_r / 2, i_theta / 2);
                fine_values[fine_idx] = (4.0 * fine_values[fine_idx] - coarse_values[coarse_idx]) / 3.0;
            }
        });

    /* Radial Indexing Section */
    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Extrapolation: Apply Extrapolation (Radial)",
        Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, fineGrid.numberSmootherCircles()},
                                                                              {fineGrid.ntheta(), fineGrid.nr()}),
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            const int fine_idx = fineGrid.index(i_r, i_theta);
            if (i_r & 1 || i_theta & 1) {
                fine_values[fine_idx] *= 4.0 / 3.0;
            }
            else {
                const int coarse_idx  = coarseGrid.index(i_r / 2, i_theta / 2);
                fine_values[fine_idx] = (4.0 * fine_values[fine_idx] - coarse_values[coarse_idx]) / 3.0;
            }
        });

    Kokkos::fence();
}

// =============================================================================
//   RHS Hierarchy
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::initRhsHierarchy(HostVector<double> rhs)
{
    Kokkos::deep_copy(levels_[0].rhs(), rhs);
    for (int level_depth = 0; level_depth < number_of_levels_ - 1; ++level_depth) {
        Level<DomainGeometry, DensityProfileCoefficients>& current_level = levels_[level_depth];
        Level<DomainGeometry, DensityProfileCoefficients>& next_level    = levels_[level_depth + 1];
        restriction(level_depth, next_level.rhs(), current_level.rhs());
    }
}

// =============================================================================
//   Convergence and Error Analysis Functions
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
std::pair<double, double> GMGPolar<DomainGeometry, DensityProfileCoefficients>::evaluateExactError(
    const PolarGrid<DefaultMemorySpace>& grid, HostConstVector<double> discrete_solution,
    HostConstVector<double> analytical_solution_host, HostVector<double> error)
{
    // Transfer the exact solution values from host to device memory for error computation.
    Kokkos::deep_copy(error, analytical_solution_host_);

    // Compute the error as the difference between the exact and numerical solutions.
    subtract(error, discrete_solution);

    // Compute the weighted L2 norm and infinity norm of the error between the numerical and exact solution.
    // The results are stored as a pair: (weighted L2 error, infinity error).
    const double weighted_euclidean_error = l2_norm(HostConstVector<double>(error)) / std::sqrt(grid.numberOfNodes());
    const double infinity_error           = infinity_norm(HostConstVector<double>(error));

    return std::make_pair(weighted_euclidean_error, infinity_error);
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::computeAnalyticalSolutionOnHost(
    const PolarGrid<Kokkos::HostSpace>& grid, HostVector<double> analytical_solution_host,
    const ExactSolution& exact_solution)
{
    assert(std::ssize(analytical_solution_host) == grid.numberOfNodes());

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Analytical Solution: Compute Values on Host (Circular)",
        Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>( // Host parallel loop
            {0, 0}, // Starting index
            {grid.numberSmootherCircles(), grid.ntheta()}), // Ending index
        [&](const int i_r, const int i_theta) {
            const int index                 = grid.index(i_r, i_theta);
            const double radius             = grid.radius(i_r);
            const double theta              = grid.theta(i_theta);
            analytical_solution_host[index] = exact_solution.exact_solution(radius, theta);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Analytical Solution: Compute Values on Host (Radial)",
        Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>( // Host parallel loop
            {0, grid.numberSmootherCircles()}, // Starting index
            {grid.ntheta(), grid.nr()}), // Ending index
        [&](const int i_theta, const int i_r) {
            const int index                 = grid.index(i_r, i_theta);
            const double radius             = grid.radius(i_r);
            const double theta              = grid.theta(i_theta);
            analytical_solution_host[index] = exact_solution.exact_solution(radius, theta);
        });

    /* Kokkos::fence() is not needed here since everything is computed on the host */
}

// =============================================================================
//   Convergence Functions
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
bool GMGPolar<DomainGeometry, DensityProfileCoefficients>::converged(double residual_norm,
                                                                     double relative_residual_norm)
{
    if (relative_tolerance_.has_value()) {
        if (!(relative_residual_norm > relative_tolerance_.value())) {
            return true;
        }
    }
    if (absolute_tolerance_.has_value()) {
        if (!(residual_norm > absolute_tolerance_.value())) {
            return true;
        }
    }
    return false;
}

// =============================================================================
//   Output and Logging Functions
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::printIterationHeader(bool is_exact_solution_provided)
{
    if (verbose_ <= 0)
        return;

    const int table_spacing = 4;
    std::cout << "------------------------------\n";
    std::cout << "---- Multigrid Iterations ----\n";
    std::cout << "------------------------------\n";
    std::cout << std::left;
    std::cout << std::setw(3 + table_spacing) << "it";
    if (absolute_tolerance_.has_value() || relative_tolerance_.has_value()) {
        std::cout << std::setw(9 + table_spacing) << "||r_k||";
        if (!FMG_)
            std::cout << std::setw(15 + table_spacing) << "||r_k||/||r_0||";
        else
            std::cout << std::setw(15 + table_spacing) << "||r_k||/||rhs||";
    }
    if (is_exact_solution_provided) {
        std::cout << std::setw(12 + table_spacing) << "||u-u_k||_l2";
        std::cout << std::setw(13 + table_spacing) << "||u-u_k||_inf";
    }
    std::cout << "\n";
    std::cout << std::right << std::setfill(' ');
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::printIterationInfo(int iteration,
                                                                              double current_residual_norm,
                                                                              double current_relative_residual_norm,
                                                                              bool is_exact_solution_provided)
{
    if (verbose_ <= 0)
        return;

    const int table_spacing = 4;
    std::cout << std::left << std::scientific << std::setprecision(2);
    std::cout << std::setw(3 + table_spacing) << iteration;
    if (absolute_tolerance_.has_value() || relative_tolerance_.has_value()) {
        std::cout << std::setw(9 + table_spacing + 2) << current_residual_norm;
        std::cout << std::setw(15 + table_spacing) << current_relative_residual_norm;
    }
    if (is_exact_solution_provided) {
        std::cout << std::setw(12 + table_spacing) << exact_errors_.back().first;
        std::cout << std::setw(13 + table_spacing) << exact_errors_.back().second;
    }
    std::cout << "\n";
    std::cout << std::right << std::defaultfloat << std::setprecision(6) << std::setfill(' ');
}
