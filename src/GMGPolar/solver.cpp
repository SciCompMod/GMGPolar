#include "../../include/GMGPolar/gmgpolar.h"

// =============================================================================
//   Full Multigrid Approximation
// =============================================================================

void IGMGPolar::fullMultigridApproximation(MultigridCycleType FMG_cycle, int FMG_iterations)
{
    // Start from the coarsest level
    int coarsest_depth    = number_of_levels_ - 1;
    Level& coarsest_level = levels_[coarsest_depth];

    // Solve directly on the coarsest level
    Kokkos::deep_copy(coarsest_level.solution(), coarsest_level.rhs());
    coarsest_level.directSolveInPlace(coarsest_level.solution()); // Direct solve on coarsest grid

    // Prolongate the solution from the coarsest level up to the finest, while applying Multigrid Cycles on each level
    for (int depth = coarsest_depth; depth > 0; --depth) {
        Level& coarse_level = levels_[depth]; // Current coarse level
        Level& fine_level   = levels_[depth - 1]; // Next finer level

        // The bi-cubic FMG interpolation is of higher order
        FMGInterpolation(coarse_level.level_depth(), fine_level.solution(), coarse_level.solution());

        // Apply some FMG iterations, except on the finest level,
        // where the interpolated solution is sufficintly accurate as an initial guess
        if (fine_level.level_depth() > 0) {
            applyMultigridIterations(fine_level, FMG_cycle, FMG_iterations);
        }
    }
}

// =============================================================================
//   Residual Handling Functions
// =============================================================================

void IGMGPolar::updateResidualNorms(Level& level, int iteration, double& initial_residual_norm,
                                    double& current_residual_norm, double& current_relative_residual_norm)
{
    level.computeResidual(level.residual(), level.rhs(), level.solution());
    if (extrapolation_ != ExtrapolationType::NONE) {
        Level& next_level = levels_[level.level_depth() + 1];
        injection(level.level_depth(), next_level.solution(), level.solution());
        next_level.computeResidual(next_level.residual(), next_level.rhs(), next_level.solution());
        extrapolatedResidual(level.level_depth(), level.residual(), next_level.residual());
    }

    current_residual_norm = residualNorm(residual_norm_type_, level, level.residual());
    residual_norms_.push_back(current_residual_norm);

    if (number_of_iterations_ == 0) {
        initial_residual_norm = !FMG_ ? current_residual_norm : residualNorm(residual_norm_type_, level, level.rhs());
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

// =============================================================================
//   Convergence and Error Analysis Functions
// =============================================================================

bool IGMGPolar::converged(double residual_norm, double relative_residual_norm)
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

void IGMGPolar::printIterationHeader(const ExactSolution* exact_solution)
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

void IGMGPolar::printIterationInfo(int iteration, double current_residual_norm, double current_relative_residual_norm,
                                   const ExactSolution* exact_solution)
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
    if (exact_solution != nullptr) {
        std::cout << std::setw(12 + table_spacing) << exact_errors_.back().first;
        std::cout << std::setw(13 + table_spacing) << exact_errors_.back().second;
    }
    std::cout << "\n";
    std::cout << std::right << std::defaultfloat << std::setprecision(6) << std::setfill(' ');
}

void IGMGPolar::solveMultigrid(double& initial_residual_norm, double& current_residual_norm,
                               double& current_relative_residual_norm)
{
    Level& level = levels_[0];

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
void IGMGPolar::solvePCG(double& initial_residual_norm, double& current_residual_norm,
                         double& current_relative_residual_norm)
{
    Level& level = levels_[0];

    // x = initial guess
    Kokkos::deep_copy(pcg_solution_, level.solution());
    // r = residual of initial guess
    Kokkos::deep_copy(level.rhs(), level.residual());

    // z = M^{-1} * r (preconditioned residual)
    if (PCG_FMG_) {
        initRhsHierarchy(level.rhs());
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
    double r_z = dot_product(ConstVector<double>(level.rhs()), ConstVector<double>(level.solution()));

    while (number_of_iterations_ < max_iterations_) {

        // A_p = A * p
        level.applySystemOperator(level.residual(), pcg_search_direction_);
        if (extrapolation_ != ExtrapolationType::NONE) {
            assert(number_of_levels_ > 1);
            Level& next_level = levels_[level.level_depth() + 1];
            injection(0, next_level.solution(), pcg_search_direction_);
            next_level.applySystemOperator(next_level.residual(), next_level.solution());
            extrapolatedResidual(0, level.residual(), next_level.residual());
        }

        // alpha = (r^T * z) / (p^T * A*p)
        double alpha =
            r_z / dot_product(ConstVector<double>(pcg_search_direction_), ConstVector<double>(level.residual()));

        // x += alpha * p
        linear_combination(pcg_solution_, 1.0, ConstVector<double>(pcg_search_direction_), alpha);

        // r -= alpha * A*p
        linear_combination(level.rhs(), 1.0, ConstVector<double>(level.residual()), -alpha);

        /* ---------------------------- */
        /* Compute convergence criteria */
        /* ---------------------------- */
        auto start_check_convergence = std::chrono::high_resolution_clock::now();

        current_residual_norm = residualNorm(residual_norm_type_, level, level.rhs());
        residual_norms_.push_back(current_residual_norm);
        current_relative_residual_norm = current_residual_norm / initial_residual_norm;

        auto end_check_convergence = std::chrono::high_resolution_clock::now();
        t_check_convergence_ += std::chrono::duration<double>(end_check_convergence - start_check_convergence).count();

        /* ---------------------------------------------- */
        /* Test solution against exact solution if given. */
        /* ---------------------------------------------- */
        LIKWID_STOP("Solver");
        auto start_check_exact_error = std::chrono::high_resolution_clock::now();

        if (exact_solution_ != nullptr)
            exact_errors_.push_back(computeExactError(level, pcg_solution_, level.residual(), *exact_solution_));

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
            initRhsHierarchy(level.rhs());
            fullMultigridApproximation(PCG_FMG_cycle_, PCG_FMG_iterations_);
        }
        else {
            // z = I^{-1} * r (no preconditioning)
            Kokkos::deep_copy(level.solution(), level.rhs());
        }
        applyMultigridIterations(level, PCG_MG_cycle_, PCG_MG_iterations_);

        // r^T * z
        double r_z_new = dot_product(ConstVector<double>(level.rhs()), ConstVector<double>(level.solution()));
        // beta = (current r^T * z) / (previous r^T * z)
        double beta = r_z_new / r_z;
        // r_z = r^T * z for next iteration
        r_z = r_z_new;

        // p *= beta
        multiply(pcg_search_direction_, beta);
        // p += z
        add(pcg_search_direction_, ConstVector<double>(level.solution()));
    }
    // level.solution() = x
    Kokkos::deep_copy(level.solution(), pcg_solution_);
}

void IGMGPolar::initRhsHierarchy(Vector<double> rhs)
{
    Kokkos::deep_copy(levels_[0].rhs(), rhs);
    for (int level_depth = 0; level_depth < number_of_levels_ - 1; ++level_depth) {
        Level& current_level = levels_[level_depth];
        Level& next_level    = levels_[level_depth + 1];
        restriction(level_depth, next_level.rhs(), current_level.rhs());
    }
}

void IGMGPolar::applyMultigridIterations(Level& level, MultigridCycleType cycle, int iterations)
{
    for (int i = 0; i < iterations; i++) {
        switch (cycle) {
        case MultigridCycleType::V_CYCLE:
            multigrid_V_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            break;
        case MultigridCycleType::W_CYCLE:
            multigrid_W_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            break;
        case MultigridCycleType::F_CYCLE:
            multigrid_F_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            break;
        default:
            std::cerr << "Error: Unknown multigrid cycle type!" << std::endl;
            break;
        }
    }
}

void IGMGPolar::applyExtrapolatedMultigridIterations(Level& level, MultigridCycleType cycle, int iterations)
{
    for (int i = 0; i < iterations; i++) {
        switch (cycle) {
        case MultigridCycleType::V_CYCLE:
            extrapolated_multigrid_V_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            break;
        case MultigridCycleType::W_CYCLE:
            extrapolated_multigrid_W_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            break;
        case MultigridCycleType::F_CYCLE:
            extrapolated_multigrid_F_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            break;
        default:
            std::cerr << "Error: Unknown multigrid cycle type!" << std::endl;
            break;
        }
    }
}
