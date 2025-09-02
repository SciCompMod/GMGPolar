#include "../../include/GMGPolar/gmgpolar.h"

#include <chrono>

// =============================================================================
//   Main Solver Routine
// =============================================================================

void GMGPolar::solve(const BoundaryConditions& boundary_conditions, const SourceTerm& source_term)
{
    LIKWID_STOP("Setup");
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
    LIKWID_START("Setup");

    LIKWID_START("Solve");
    auto start_solve = std::chrono::high_resolution_clock::now();

    // Clear solve-phase timings
    resetSolvePhaseTimings();

    /* ---------------------------- */
    /* Initialize starting solution */
    /* ---------------------------- */
    auto start_initial_approximation = std::chrono::high_resolution_clock::now();

    initializeSolution(boundary_conditions);

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
        if (exact_solution_ != nullptr)
            evaluateExactError(level, *exact_solution_);

        /* ---------------------------- */
        /* Compute convergence criteria */
        /* ---------------------------- */
        auto start_check_convergence = std::chrono::high_resolution_clock::now();
        if (absolute_tolerance_.has_value() || relative_tolerance_.has_value()) {
            updateResidualNorms(level, number_of_iterations_, initial_residual_norm, current_residual_norm,
                                current_relative_residual_norm);
        }

        printIterationInfo(number_of_iterations_, current_residual_norm, current_relative_residual_norm,
                           exact_solution_);

        if (converged(current_residual_norm, current_relative_residual_norm))
            break;

        auto end_check_convergence = std::chrono::high_resolution_clock::now();
        t_check_convergence_ += std::chrono::duration<double>(end_check_convergence - start_check_convergence).count();

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

// =============================================================================
//   Solution Initialization
// =============================================================================

void GMGPolar::initializeSolution(const BoundaryConditions& boundary_conditions)
{
    if (!FMG_) {
        int start_level_depth = 0;
        Level& level          = levels_[start_level_depth];
        assign(level.solution(), 0.0); // Assign zero initial guess if not using FMG

        /* Consider setting the boundary conditions u_D and u_D_Interior if DirBC_Interior to the initial solution */
        bool use_boundary_condition = false;
        if (use_boundary_condition) {
            const auto& grid            = level.grid();
            const auto& sin_theta_cache = level.levelCache().sin_theta();
            const auto& cos_theta_cache = level.levelCache().cos_theta();

            const int i_r_inner  = 0;
            const int i_r_outer  = grid.nr() - 1;
            const double r_inner = grid.radius(i_r_inner);
            const double r_outer = grid.radius(i_r_outer);

            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                const double theta     = grid.theta(i_theta);
                const double sin_theta = sin_theta_cache[i_theta];
                const double cos_theta = cos_theta_cache[i_theta];
                if (DirBC_Interior_) { // Apply interior Dirichlet BC if enabled.
                    const int index         = grid.index(i_r_inner, i_theta);
                    level.solution()[index] = boundary_conditions.u_D_Interior(r_inner, theta, sin_theta, cos_theta);
                }
                // Always apply outer boundary condition.
                const int index         = grid.index(i_r_outer, i_theta);
                level.solution()[index] = boundary_conditions.u_D(r_outer, theta, sin_theta, cos_theta);
            }
        }
    }
    else {
        // Start from the coarsest level
        int FMG_start_level_depth = number_of_levels_ - 1;
        Level& FMG_level          = levels_[FMG_start_level_depth];

        // Solve directly on the coarsest level
        FMG_level.solution() = FMG_level.rhs();
        FMG_level.directSolveInPlace(FMG_level.solution()); // Direct solve on coarsest grid

        // Prolongate the solution from the coarsest level up to the finest, while applying Multigrid Cycles on each level
        for (int current_level = FMG_start_level_depth - 1; current_level > 0; --current_level) {
            Level& FMG_level      = levels_[current_level]; // The current level
            Level& next_FMG_level = levels_[current_level - 1]; // The finer level

            // The bi-cubic FMG interpolation is of higher order
            FMGInterpolation(current_level, next_FMG_level.solution(), FMG_level.solution());

            // Apply some FMG iterations
            for (int i = 0; i < FMG_iterations_; i++) {
                if (current_level - 1 == 0 && (extrapolation_ != ExtrapolationType::NONE)) {
                    switch (FMG_cycle_) {
                    case MultigridCycleType::V_CYCLE:
                        implicitlyExtrapolatedMultigrid_V_Cycle(current_level - 1, next_FMG_level.solution(),
                                                                next_FMG_level.rhs(), next_FMG_level.residual());
                        break;

                    case MultigridCycleType::W_CYCLE:
                        implicitlyExtrapolatedMultigrid_W_Cycle(current_level - 1, next_FMG_level.solution(),
                                                                next_FMG_level.rhs(), next_FMG_level.residual());
                        break;

                    case MultigridCycleType::F_CYCLE:
                        implicitlyExtrapolatedMultigrid_F_Cycle(current_level - 1, next_FMG_level.solution(),
                                                                next_FMG_level.rhs(), next_FMG_level.residual());
                        break;

                    default:
                        std::cerr << "Error: Unknown multigrid cycle type!" << std::endl;
                        throw std::runtime_error("Invalid multigrid cycle type encountered.");
                        break;
                    }
                }
                else {
                    switch (FMG_cycle_) {
                    case MultigridCycleType::V_CYCLE:
                        multigrid_V_Cycle(current_level - 1, next_FMG_level.solution(), next_FMG_level.rhs(),
                                          next_FMG_level.residual());
                        break;

                    case MultigridCycleType::W_CYCLE:
                        multigrid_W_Cycle(current_level - 1, next_FMG_level.solution(), next_FMG_level.rhs(),
                                          next_FMG_level.residual());
                        break;

                    case MultigridCycleType::F_CYCLE:
                        multigrid_F_Cycle(current_level - 1, next_FMG_level.solution(), next_FMG_level.rhs(),
                                          next_FMG_level.residual());
                        break;

                    default:
                        std::cerr << "Error: Unknown multigrid cycle type!" << std::endl;
                        throw std::runtime_error("Invalid multigrid cycle type encountered.");
                        break;
                    }
                }
            }
        }
    }
}

// =============================================================================
//   Residual Handling Functions
// =============================================================================

void GMGPolar::updateResidualNorms(Level& level, int iteration, double& initial_residual_norm,
                                   double& current_residual_norm, double& current_relative_residual_norm)
{
    level.computeResidual(level.residual(), level.rhs(), level.solution());

    if (extrapolation_ != ExtrapolationType::NONE) {
        Level& next_level = levels_[level.level_depth() + 1];
        injection(level.level_depth(), next_level.solution(), level.solution());
        next_level.computeResidual(next_level.residual(), next_level.rhs(), next_level.solution());
        extrapolatedResidual(level.level_depth(), level.residual(), next_level.residual());
    }

    switch (residual_norm_type_) {
    case ResidualNormType::EUCLIDEAN:
        current_residual_norm = l2_norm(level.residual());
        break;
    case ResidualNormType::WEIGHTED_EUCLIDEAN:
        current_residual_norm = l2_norm(level.residual()) / std::sqrt(level.grid().numberOfNodes());
        break;
    case ResidualNormType::INFINITY_NORM:
        current_residual_norm = infinity_norm(level.residual());
        break;
    default:
        throw std::invalid_argument("Unknown ResidualNormType");
    }
    residual_norms_.push_back(current_residual_norm);

    if (number_of_iterations_ == 0) {
        if (!FMG_) {
            initial_residual_norm = current_residual_norm;
        }
        else {
            switch (residual_norm_type_) {
            case ResidualNormType::EUCLIDEAN:
                initial_residual_norm = l2_norm(level.rhs());
                break;
            case ResidualNormType::WEIGHTED_EUCLIDEAN:
                initial_residual_norm = (l2_norm(level.rhs()) / std::sqrt(level.grid().numberOfNodes()));
                break;
            case ResidualNormType::INFINITY_NORM:
                initial_residual_norm = infinity_norm(level.rhs());
                break;
            default:
                throw std::invalid_argument("Unknown ResidualNormType");
            }
        }
    }

    current_relative_residual_norm = current_residual_norm / initial_residual_norm;

    if (number_of_iterations_ > 0) {
        const double current_residual_reduction_factor =
            residual_norms_[number_of_iterations_] / residual_norms_[number_of_iterations_ - 1];

        const double convergence_factor = 0.7;
        if (current_residual_reduction_factor > convergence_factor && extrapolation_ == ExtrapolationType::COMBINED &&
            full_grid_smoothing_) {
            full_grid_smoothing_ = false;
            if (verbose_ > 0)
                std::cout << "Switching from full grid smoothing to standard extrapolated smoothing.\n";
        }
    }
}

void GMGPolar::extrapolatedResidual(const int current_level, Vector<double>& residual,
                                    const Vector<double>& residual_next_level)
{
    const PolarGrid& fineGrid   = levels_[current_level].grid();
    const PolarGrid& coarseGrid = levels_[current_level + 1].grid();

    assert(residual.size() == fineGrid.numberOfNodes());
    assert(residual_next_level.size() == coarseGrid.numberOfNodes());

#pragma omp parallel num_threads(threads_per_level_[current_level])
    {
/* Circluar Indexing Section */
/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_r = 0; i_r < fineGrid.numberSmootherCircles(); i_r++) {
            int i_r_coarse = i_r / 2;
            for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++) {
                int i_theta_coarse = i_theta / 2;

                if (i_r & 1 || i_theta & 1) {
                    residual[fineGrid.index(i_r, i_theta)] *= 4.0 / 3.0;
                }
                else {
                    int fine_idx       = fineGrid.index(i_r, i_theta);
                    int coarse_idx     = coarseGrid.index(i_r_coarse, i_theta_coarse);
                    residual[fine_idx] = (4.0 * residual[fine_idx] - residual_next_level[coarse_idx]) / 3.0;
                }
            }
        }

/* Radial Indexing Section */
/* For loop matches radial access pattern */
#pragma omp for nowait
        for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++) {
            int i_theta_coarse = i_theta / 2;
            for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++) {
                int i_r_coarse = i_r / 2;

                if (i_r & 1 || i_theta & 1) {
                    residual[fineGrid.index(i_r, i_theta)] *= 4.0 / 3.0;
                }
                else {
                    int fine_idx       = fineGrid.index(i_r, i_theta);
                    int coarse_idx     = coarseGrid.index(i_r_coarse, i_theta_coarse);
                    residual[fine_idx] = (4.0 * residual[fine_idx] - residual_next_level[coarse_idx]) / 3.0;
                }
            }
        }
    }
}

// =============================================================================
//   Convergence and Error Analysis Functions
// =============================================================================

bool GMGPolar::converged(const double& residual_norm, const double& relative_residual_norm)
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

void GMGPolar::evaluateExactError(Level& level, const ExactSolution& exact_solution)
{
    LIKWID_STOP("Solver");
    auto start = std::chrono::high_resolution_clock::now();
    // Compute the weighted L2 norm and infinity norm of the error between the numerical and exact solution.
    // The results are stored as a pair: (weighted L2 error, infinity error).
    exact_errors_.push_back(computeExactError(level, level.solution(), level.residual(), exact_solution));
    auto end = std::chrono::high_resolution_clock::now();
    t_check_exact_error_ += std::chrono::duration<double>(end - start).count();
    LIKWID_START("Solver");
}

std::pair<double, double> GMGPolar::computeExactError(Level& level, const Vector<double>& solution,
                                                      Vector<double>& error, const ExactSolution& exact_solution)
{
    const PolarGrid& grid        = level.grid();
    const LevelCache& levelCache = level.levelCache();
    const auto& sin_theta_cache  = levelCache.sin_theta();
    const auto& cos_theta_cache  = levelCache.cos_theta();

    assert(solution.size() == error.size());
    assert(solution.size() == grid.numberOfNodes());

#pragma omp parallel num_threads(threads_per_level_[level.level_depth()])
    {
#pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                double theta     = grid.theta(i_theta);
                double sin_theta = sin_theta_cache[i_theta];
                double cos_theta = cos_theta_cache[i_theta];
                error[grid.index(i_r, i_theta)] =
                    exact_solution.exact_solution(r, theta, sin_theta, cos_theta) - solution[grid.index(i_r, i_theta)];
            }
        }
#pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            double theta     = grid.theta(i_theta);
            double sin_theta = sin_theta_cache[i_theta];
            double cos_theta = cos_theta_cache[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                double r = grid.radius(i_r);
                error[grid.index(i_r, i_theta)] =
                    exact_solution.exact_solution(r, theta, sin_theta, cos_theta) - solution[grid.index(i_r, i_theta)];
            }
        }
    }

    double weighted_euclidean_error = l2_norm(error) / std::sqrt(grid.numberOfNodes());
    double infinity_error           = infinity_norm(error);

    return std::make_pair(weighted_euclidean_error, infinity_error);
}

// =============================================================================
//   Output and Logging Functions
// =============================================================================

void GMGPolar::printIterationHeader(const ExactSolution* exact_solution)
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
    if (exact_solution != nullptr) {
        std::cout << std::setw(12 + table_spacing) << "||u-u_k||_l2";
        std::cout << std::setw(13 + table_spacing) << "||u-u_k||_inf";
    }
    std::cout << "\n";
    std::cout << std::right << std::setfill(' ');
}

void GMGPolar::printIterationInfo(int iteration, double current_residual_norm, double current_relative_residual_norm,
                                  const ExactSolution* exact_solution)
{
    if (verbose_ <= 0)
        return;

    const int table_spacing = 4;
    std::cout << std::left << std::scientific << std::setprecision(2);
    std::cout << std::setw(3 + table_spacing) << number_of_iterations_;
    if (absolute_tolerance_.has_value() || relative_tolerance_.has_value()) {
        std::cout << std::setw(9 + table_spacing + 2) << current_residual_norm;
        std::cout << std::setw(15 + table_spacing) << current_relative_residual_norm;
    }
    if (exact_solution != nullptr) {
        std::cout << std::setw(12 + table_spacing) << exact_errors_.back().first;
        std::cout << std::setw(13 + table_spacing) << exact_errors_.back().second;
    }
    std::cout << "\n";
    std::cout << std::right << std::defaultfloat << std::setprecision(6) << std::setfill(' ');
}
