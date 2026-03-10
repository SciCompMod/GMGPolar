// =============================================================================
//   Solution Access
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
Vector<double> GMGPolar<DomainGeometry, DensityProfileCoefficients>::solution()
{
    int level_depth = 0;
    return levels_[level_depth].solution();
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
ConstVector<double> GMGPolar<DomainGeometry, DensityProfileCoefficients>::solution() const
{
    int level_depth = 0;
    return levels_[level_depth].solution();
}

// =============================================================================
//   Solution Initialization
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::initializeSolution()
{
    if (!FMG_) {
        int start_level_depth        = 0;
        Level<DomainGeometry>& level = levels_[start_level_depth];
        assign(level.solution(), 0.0); // Assign zero initial guess if not using FMG
    }
    else {
        // Start from the coarsest level
        int coarsest_depth                    = number_of_levels_ - 1;
        Level<DomainGeometry>& coarsest_level = levels_[coarsest_depth];

        // Solve directly on the coarsest level
        Kokkos::deep_copy(coarsest_level.solution(), coarsest_level.rhs());
        coarsest_level.directSolveInPlace(coarsest_level.solution()); // Direct solve on coarsest grid

        // Prolongate the solution from the coarsest level up to the finest, while applying Multigrid Cycles on each level
        for (int depth = coarsest_depth; depth > 0; --depth) {
            Level<DomainGeometry>& coarse_level = levels_[depth];     // Current coarse level
            Level<DomainGeometry>& fine_level   = levels_[depth - 1]; // Next finer level

            // The bi-cubic FMG interpolation is of higher order
            FMGInterpolation(coarse_level.level_depth(), fine_level.solution(), coarse_level.solution());

            // Apply some FMG iterations
            for (int i = 0; i < FMG_iterations_; i++) {
                if (fine_level.level_depth() == 0 && (extrapolation_ != ExtrapolationType::NONE)) {
                    switch (FMG_cycle_) {
                    case MultigridCycleType::V_CYCLE:
                        extrapolated_multigrid_V_Cycle(fine_level.level_depth(), fine_level.solution(),
                                                       fine_level.rhs(), fine_level.residual());
                        break;

                    case MultigridCycleType::W_CYCLE:
                        extrapolated_multigrid_W_Cycle(fine_level.level_depth(), fine_level.solution(),
                                                       fine_level.rhs(), fine_level.residual());
                        break;

                    case MultigridCycleType::F_CYCLE:
                        extrapolated_multigrid_F_Cycle(fine_level.level_depth(), fine_level.solution(),
                                                       fine_level.rhs(), fine_level.residual());
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
                        multigrid_V_Cycle(fine_level.level_depth(), fine_level.solution(), fine_level.rhs(),
                                          fine_level.residual());
                        break;

                    case MultigridCycleType::W_CYCLE:
                        multigrid_W_Cycle(fine_level.level_depth(), fine_level.solution(), fine_level.rhs(),
                                          fine_level.residual());
                        break;

                    case MultigridCycleType::F_CYCLE:
                        multigrid_F_Cycle(fine_level.level_depth(), fine_level.solution(), fine_level.rhs(),
                                          fine_level.residual());
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

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
double GMGPolar<DomainGeometry, DensityProfileCoefficients>::residualNorm(const ResidualNormType& norm_type,
                                                                          const Level<DomainGeometry>& level,
                                                                          ConstVector<double> residual) const
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
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::updateResidualNorms(
    Level<DomainGeometry>& level, int iteration, double& initial_residual_norm, double& current_residual_norm,
    double& current_relative_residual_norm)
{
    level.computeResidual(level.residual(), level.rhs(), level.solution());
    if (extrapolation_ != ExtrapolationType::NONE) {
        Level<DomainGeometry>& next_level = levels_[level.level_depth() + 1];
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

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::extrapolatedResidual(int current_level,
                                                                                Vector<double> residual,
                                                                                ConstVector<double> residual_next_level)
{
    const PolarGrid& fineGrid   = levels_[current_level].grid();
    const PolarGrid& coarseGrid = levels_[current_level + 1].grid();

    assert(std::ssize(residual) == fineGrid.numberOfNodes());
    assert(std::ssize(residual_next_level) == coarseGrid.numberOfNodes());

#pragma omp parallel num_threads(max_omp_threads_)
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

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::evaluateExactError(Level<DomainGeometry>& level,
                                                                              const ExactSolution& exact_solution)
{
    // Compute the weighted L2 norm and infinity norm of the error between the numerical and exact solution.
    // The results are stored as a pair: (weighted L2 error, infinity error).
    exact_errors_.push_back(computeExactError(level, level.solution(), level.residual(), exact_solution));
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
std::pair<double, double> GMGPolar<DomainGeometry, DensityProfileCoefficients>::computeExactError(
    Level<DomainGeometry>& level, ConstVector<double> solution, Vector<double> error,
    const ExactSolution& exact_solution)
{
    const PolarGrid& grid                          = level.grid();
    const LevelCache<DomainGeometry>& levelCache   = level.levelCache();

    assert(solution.size() == error.size());
    assert(std::ssize(solution) == grid.numberOfNodes());

#pragma omp parallel num_threads(max_omp_threads_)
    {
#pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                double theta = grid.theta(i_theta);
                error[grid.index(i_r, i_theta)] =
                    exact_solution.exact_solution(r, theta) - solution[grid.index(i_r, i_theta)];
            }
        }
#pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            double theta = grid.theta(i_theta);
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                double r = grid.radius(i_r);
                error[grid.index(i_r, i_theta)] =
                    exact_solution.exact_solution(r, theta) - solution[grid.index(i_r, i_theta)];
            }
        }
    }
    ConstVector<double> c_error     = error;
    double weighted_euclidean_error = l2_norm(c_error) / std::sqrt(grid.numberOfNodes());
    double infinity_error           = infinity_norm(c_error);

    return std::make_pair(weighted_euclidean_error, infinity_error);
}

// =============================================================================
//   Full Multigrid Approximation
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::fullMultigridApproximation(MultigridCycleType FMG_cycle,
                                                                                      int FMG_iterations)
{
    // Start from the coarsest level
    int coarsest_depth                    = number_of_levels_ - 1;
    Level<DomainGeometry>& coarsest_level = levels_[coarsest_depth];

    // Solve directly on the coarsest level
    Kokkos::deep_copy(coarsest_level.solution(), coarsest_level.rhs());
    coarsest_level.directSolveInPlace(coarsest_level.solution());

    // Prolongate the solution from the coarsest level up to the finest, while applying Multigrid Cycles on each level
    for (int depth = coarsest_depth; depth > 0; --depth) {
        Level<DomainGeometry>& coarse_level = levels_[depth];     // Current coarse level
        Level<DomainGeometry>& fine_level   = levels_[depth - 1]; // Next finer level

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
//   RHS Hierarchy
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::initRhsHierarchy(Vector<double> rhs)
{
    Kokkos::deep_copy(levels_[0].rhs(), rhs);
    for (int level_depth = 0; level_depth < number_of_levels_ - 1; ++level_depth) {
        Level<DomainGeometry>& current_level = levels_[level_depth];
        Level<DomainGeometry>& next_level    = levels_[level_depth + 1];
        restriction(level_depth, next_level.rhs(), current_level.rhs());
    }
}

// =============================================================================
//   Apply Multigrid Iterations
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::applyMultigridIterations(Level<DomainGeometry>& level,
                                                                                    MultigridCycleType cycle,
                                                                                    int iterations)
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

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::applyExtrapolatedMultigridIterations(
    Level<DomainGeometry>& level, MultigridCycleType cycle, int iterations)
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

// =============================================================================
//   Solver Loops
// =============================================================================

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::solveMultigrid(double& initial_residual_norm,
                                                                          double& current_residual_norm,
                                                                          double& current_relative_residual_norm)
{
    Level<DomainGeometry>& level = levels_[0];

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
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::solvePCG(double& initial_residual_norm,
                                                                    double& current_residual_norm,
                                                                    double& current_relative_residual_norm)
{
    Level<DomainGeometry>& level = levels_[0];

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
            Level<DomainGeometry>& next_level = levels_[level.level_depth() + 1];
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

// =============================================================================
//   Convergence and Output (IGMGPolar methods)
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

// =============================================================================
//   (Legacy: do_solve replaced by solve<BC> in solver.h — kept as a stub)
// =============================================================================
// REMOVED: do_solve(BoundaryConditionRef, SourceTerm) — see solver.h for the
// templated solve<BoundaryConditions>(boundary_conditions, source_term).
#if 0
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::do_solve(const BoundaryConditionRef& boundary_conditions,
                                                                    const SourceTerm& source_term)
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
        Level<DomainGeometry>& current_level = levels_[level_depth];
        // Inject rhs if there is a next level
        if (level_depth + 1 < initial_rhs_f_levels) {
            Level<DomainGeometry>& next_level = levels_[level_depth + 1];
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

    initializeSolution();

    auto end_initial_approximation = std::chrono::high_resolution_clock::now();
    t_solve_initial_approximation_ =
        std::chrono::duration<double>(end_initial_approximation - start_initial_approximation).count();

    // These times are included in the initial approximation and don't count towards the multigrid cyclces.
    resetAvgMultigridCycleTimings();

    /* --------------------------------------- */
    /* Start Solver at finest level (depth 0)  */
    /* --------------------------------------- */
    Level<DomainGeometry>& level = levels_[0];

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
                extrapolated_multigrid_V_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            }
            break;
        case MultigridCycleType::W_CYCLE:
            if (extrapolation_ == ExtrapolationType::NONE) {
                multigrid_W_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            }
            else {
                extrapolated_multigrid_W_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            }
            break;
        case MultigridCycleType::F_CYCLE:
            if (extrapolation_ == ExtrapolationType::NONE) {
                multigrid_F_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
            }
            else {
                extrapolated_multigrid_F_Cycle(level.level_depth(), level.solution(), level.rhs(), level.residual());
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
#endif
