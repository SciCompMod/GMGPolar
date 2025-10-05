
#include "../../include/GMGPolar/gmgpolar.h"

#include <random>

void GMGPolar::solvePCG(const BoundaryConditions& boundary_conditions, const SourceTerm& source_term,
                        int PCG_FMG_iterations, MultigridCycleType PCG_FMG_cycle, ExtrapolationType PCG_extrapolation)
{

    std::cout << "------------------------------\n";
    std::cout << "---- Preconditioner (CG) -----\n";
    std::cout << "------------------------------\n";

    switch (PCG_extrapolation) {
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

    std::cout << "Full-Multigrid: " << PCG_FMG_iterations << "x ";
    switch (PCG_FMG_cycle) {
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

    /* --------------------------------------- */
    /* Start Solver at finest level (depth 0)  */
    /* --------------------------------------- */
    Level& level_0 = levels_[0];
    Level& level_1 = levels_[1];

    /* Define Preconditioner */
    GMGPolar preconditioner(level_0.grid(), domain_geometry_, density_profile_coefficients_);
    preconditioner.verbose(0); // Enable/disable verbose output
    preconditioner.paraview(0); // Enable/disable ParaView output
    preconditioner.maxOpenMPThreads(max_omp_threads_); // Maximum OpenMP threads to use
    preconditioner.threadReductionFactor(thread_reduction_factor_); // Reduce threads on coarser grids
    preconditioner.DirBC_Interior(DirBC_Interior_); // Interior boundary conditions: Dirichlet, Across-the-origin,
    preconditioner.stencilDistributionMethod(stencil_distribution_method_); // Stencil distribution strategy: Take, Give
    preconditioner.cacheDensityProfileCoefficients(
        cache_density_profile_coefficients_); // Cache density profile coefficients: alpha, beta
    preconditioner.cacheDomainGeometry(cache_domain_geometry_); // Cache domain geometry data: arr, att, art, detDF
    preconditioner.extrapolation(PCG_extrapolation); // Enable/disable extrapolation
    preconditioner.maxLevels(max_levels_); // Max multigrid levels (-1 = use deepest possible)
    preconditioner.preSmoothingSteps(pre_smoothing_steps_); // Smoothing before coarse-grid correction
    preconditioner.postSmoothingSteps(post_smoothing_steps_); // Smoothing after coarse-grid correction
    preconditioner.multigridCycle(multigrid_cycle_); // Multigrid cycle type
    preconditioner.FMG(true);
    preconditioner.FMG_iterations(PCG_FMG_iterations); // FMG iteration count
    preconditioner.FMG_cycle(PCG_FMG_cycle); // FMG cycle type
    preconditioner.maxIterations(0); // Max number of iterations
    preconditioner.residualNormType(residual_norm_type_); // Residual norm type (L2, weighted-L2, L∞)
    preconditioner.absoluteTolerance(absolute_tolerance_); // Absolute residual tolerance
    preconditioner.relativeTolerance(relative_tolerance_); // Relative residual tolerance
    preconditioner.setup(); // (allocates internal data, prepares operators, etc.)

    number_of_iterations_                 = 0;
    double initial_residual_norm          = 1.0;
    double current_residual_norm          = 1.0;
    double current_relative_residual_norm = 1.0;

    const int n_0 = level_0.grid().numberOfNodes();
    const int n_1 = level_1.grid().numberOfNodes();
    Vector<double> z_0(n_0);
    Vector<double> p_0(n_0);
    Vector<double> A_p_0(n_0);
    Vector<double> z_1(n_1);
    Vector<double> p_1(n_1);
    Vector<double> A_p_1(n_1);
    Vector<double> error(n_0);

    LIKWID_START("Solve");
    auto start_solve = std::chrono::high_resolution_clock::now();

    // Clear solve-phase timings
    resetSolvePhaseTimings();

    /* ---------------------------- */
    /* Initialize starting solution */
    /* ---------------------------- */
    auto start_initial_approximation = std::chrono::high_resolution_clock::now();

    initializeSolution();

    auto end_initial_approximation = std::chrono::high_resolution_clock::now();
    t_solve_initial_approximation_ =
        std::chrono::duration<double>(end_initial_approximation - start_initial_approximation).count();

    // These times are included in the initial approximation and don't count towards the multigrid cyclces.
    resetAvgMultigridCycleTimings();

    printIterationHeader(exact_solution_);

    // Compute Residual for the initial guess
    level_0.computeResidual(level_0.residual(), level_0.rhs(), level_0.solution());
    if (extrapolation_ != ExtrapolationType::NONE) {
        injection(0, level_1.solution(), level_0.solution());
        level_1.computeResidual(level_1.residual(), level_1.rhs(), level_1.solution());
        extrapolatedResidual(0, level_0.residual(), level_1.residual());
    }

    current_residual_norm = residualNorm(residual_norm_type_, level_0, level_0.residual());
    residual_norms_.push_back(current_residual_norm);
    if (number_of_iterations_ == 0) {
        initial_residual_norm =
            !FMG_ ? current_residual_norm : residualNorm(residual_norm_type_, level_0, level_0.rhs());
    }
    current_relative_residual_norm = current_residual_norm / initial_residual_norm;

    // --- Check exact error, excluded from timings
    auto start_check_exact_error = std::chrono::high_resolution_clock::now();
    if (exact_solution_ != nullptr)
        exact_errors_.push_back(computeExactError(level_0, level_0.solution(), error, *exact_solution_));
    auto end_check_exact_error = std::chrono::high_resolution_clock::now();
    t_check_exact_error_ += std::chrono::duration<double>(end_check_exact_error - start_check_exact_error).count();

    LIKWID_START("Solver");

    printIterationInfo(number_of_iterations_, current_residual_norm, current_relative_residual_norm, exact_solution_);

    /* ---------------------------- */
    /* Start solver timing (CG part)*/
    /* ---------------------------- */
    auto start_cg_solve = std::chrono::high_resolution_clock::now();

    preconditioner.solve(level_0.residual());
    z_0 = preconditioner.solution();
    p_0 = z_0;

    double r_z = dot_product(level_0.residual(), z_0);

    while (number_of_iterations_ < max_iterations_) {
        // A_p = A * p
        level_0.applySystemOperator(A_p_0, p_0);
        if (extrapolation_ != ExtrapolationType::NONE) {
            injection(0, p_1, p_0);
            level_1.applySystemOperator(A_p_1, p_1);
            extrapolatedResidual(0, A_p_0, A_p_1);
        }

        double alpha = r_z / dot_product(p_0, A_p_0);

        linear_combination(level_0.solution(), 1.0, p_0, alpha);
        linear_combination(level_0.residual(), 1.0, A_p_0, -alpha);

        current_residual_norm = residualNorm(residual_norm_type_, level_0, level_0.residual());
        residual_norms_.push_back(current_residual_norm);
        current_relative_residual_norm = current_residual_norm / initial_residual_norm;

        // --- Check exact error, excluded from timings
        auto start_check_exact_error = std::chrono::high_resolution_clock::now();
        if (exact_solution_ != nullptr)
            exact_errors_.push_back(computeExactError(level_0, level_0.solution(), error, *exact_solution_));
        auto end_check_exact_error = std::chrono::high_resolution_clock::now();
        t_check_exact_error_ += std::chrono::duration<double>(end_check_exact_error - start_check_exact_error).count();

        number_of_iterations_++;

        printIterationInfo(number_of_iterations_, current_residual_norm, current_relative_residual_norm,
                           exact_solution_);

        if (converged(current_residual_norm, current_relative_residual_norm))
            break;

        preconditioner.solve(level_0.residual());
        z_0 = preconditioner.solution();

        double r_z_new = dot_product(level_0.residual(), z_0);
        double beta    = r_z_new / r_z;

        r_z = r_z_new;

        multiply(p_0, beta);
        add(p_0, z_0);
    }

    auto end_solve    = std::chrono::high_resolution_clock::now();
    auto end_cg_solve = std::chrono::high_resolution_clock::now();
    t_solve_total_    = std::chrono::duration<double>(end_solve - start_solve).count() - t_check_exact_error_;
    double t_solve_cg = std::chrono::duration<double>(end_cg_solve - start_cg_solve).count() - t_check_exact_error_;

    LIKWID_STOP("Solve");

    /* ---------------------- */
    /* Print Timing Summary   */
    /* ---------------------- */
    std::cout << "------------------\n"
              << "Timing Information\n"
              << "------------------\n"
              << "Solve Time: " << t_solve_total_ << " seconds\n"
              << "    Initial Approximation: " << t_solve_initial_approximation_ << " seconds\n"
              << "    CG Iteration: " << t_solve_cg << " seconds\n";

    if (number_of_iterations_ > 0) {
        std::cout << "Average CG Iteration: " << (t_solve_cg / number_of_iterations_) << " seconds\n";
    }
}

void GMGPolar::solve(const Vector<double>& rhs_f)
{
    auto start_setup_rhs = std::chrono::high_resolution_clock::now();

    /* ------------------------------------- */
    /* Build rhs_f on Level 0 (finest Level) */
    /* ------------------------------------- */
    levels_[0].rhs() = rhs_f;

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
            restriction(level_depth, next_level.rhs(), current_level.rhs());
        }
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

    initializeSolution();

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

    while (number_of_iterations_ < max_iterations_) {
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

        printIterationInfo(number_of_iterations_, current_residual_norm, current_relative_residual_norm,
                           exact_solution_);

        if (converged(current_residual_norm, current_relative_residual_norm))
            break;

        /* ----------------------- */
        /* Perform Multigrid Cycle */
        /* ----------------------- */
        auto start_solve_multigrid_iterations = std::chrono::high_resolution_clock::now();

        switch (multigrid_cycle_) {
        case MultigridCycleType::V_CYCLE:
            if (extrapolation_ == ExtrapolationType::NONE) {
                multigrid_V_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            }
            else {
                implicitlyExtrapolatedMultigrid_V_Cycle(level.level_depth(), level.solution(), level.rhs(),
                                                        level.residual());
            }
            break;
        case MultigridCycleType::W_CYCLE:
            if (extrapolation_ == ExtrapolationType::NONE) {
                multigrid_W_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            }
            else {
                implicitlyExtrapolatedMultigrid_W_Cycle(level.level_depth(), level.solution(), level.rhs(),
                                                        level.residual());
            }
            break;
        case MultigridCycleType::F_CYCLE:
            if (extrapolation_ == ExtrapolationType::NONE) {
                multigrid_F_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            }
            else {
                implicitlyExtrapolatedMultigrid_F_Cycle(level.level_depth(), level.solution(), level.rhs(),
                                                        level.residual());
            }
            break;
        default:
            throw std::invalid_argument("Unknown MultigridCycleType");
        }
        number_of_iterations_++;

        auto end_solve_multigrid_iterations = std::chrono::high_resolution_clock::now();
        t_solve_multigrid_iterations_ +=
            std::chrono::duration<double>(end_solve_multigrid_iterations - start_solve_multigrid_iterations).count();
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
