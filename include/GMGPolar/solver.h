// =============================================================================
//   Main Solver Routine
// =============================================================================
template <concepts::BoundaryConditions BoundaryConditions>
void IGMGPolar::solve(const BoundaryConditions& boundary_conditions, const SourceTerm& source_term)
{
    auto start_setup_rhs = std::chrono::high_resolution_clock::now();

    /* ------------------------------------- */
    /* Build rhs_f on Level 0 (finest Level) */
    /* ------------------------------------- */
    build_rhs_f(levels_[0], levels_[0].rhs(), boundary_conditions, source_term);

    /* ---------------- */
    /* Discretize rhs_f */
    /* ---------------- */
    int initial_rhs_f_levels = FMG_ ? number_of_levels_ : (extrapolation_ == ExtrapolationType::NONE ? 1 : 2);
    // Loop through the levels, injecting and discretizing rhs
    for (int level_depth = 0; level_depth < initial_rhs_f_levels; ++level_depth) {
        Level& current_level = levels_[level_depth];
        // Inject rhs if there is a next level
        if (level_depth + 1 < initial_rhs_f_levels) {
            Level& next_level = levels_[level_depth + 1];
            injection(level_depth, next_level.rhs(), current_level.rhs());
        }
        // Discretize the rhs for the current level
        discretize_rhs_f(current_level, current_level.rhs());
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
    Level& level = levels_[0];

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

    if (exact_solution_ != nullptr)
        evaluateExactError(level, *exact_solution_);

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
            solvePCG(initial_residual_norm, current_residual_norm, current_relative_residual_norm);
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
        if (exact_solution_ != nullptr) {
            computeExactError(level, level.solution(), level.residual(), *exact_solution_);
            writeToVTK("output_error", level, level.residual());
        }
    }
}
