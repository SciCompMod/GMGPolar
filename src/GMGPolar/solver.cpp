#include "../../include/GMGPolar/gmgpolar.h"

#include <chrono>

// =============================================================================
//   Solution Initialization
// =============================================================================

void IGMGPolar::initializeSolution()
{
    if (!FMG_) {
        int start_level_depth = 0;
        Level& level          = levels_[start_level_depth];
        assign(level.solution(), 0.0); // Assign zero initial guess if not using FMG
    }
    else {
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

            // Apply some FMG iterations
            for (int i = 0; i < FMG_iterations_; i++) {
                if (fine_level.level_depth() == 0 && (extrapolation_ != ExtrapolationType::NONE)) {
                    switch (FMG_cycle_) {
                    case MultigridCycleType::V_CYCLE:
                        implicitlyExtrapolatedMultigrid_V_Cycle(fine_level.level_depth(), fine_level.solution(),
                                                                fine_level.rhs(), fine_level.residual());
                        break;

                    case MultigridCycleType::W_CYCLE:
                        implicitlyExtrapolatedMultigrid_W_Cycle(fine_level.level_depth(), fine_level.solution(),
                                                                fine_level.rhs(), fine_level.residual());
                        break;

                    case MultigridCycleType::F_CYCLE:
                        implicitlyExtrapolatedMultigrid_F_Cycle(fine_level.level_depth(), fine_level.solution(),
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

double IGMGPolar::residualNorm(const ResidualNormType& norm_type, const Level& level,
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

void IGMGPolar::extrapolatedResidual(const int current_level, Vector<double> residual,
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

void IGMGPolar::evaluateExactError(Level& level, const ExactSolution& exact_solution)
{
    // Compute the weighted L2 norm and infinity norm of the error between the numerical and exact solution.
    // The results are stored as a pair: (weighted L2 error, infinity error).
    exact_errors_.push_back(computeExactError(level, level.solution(), level.residual(), exact_solution));
}

std::pair<double, double> IGMGPolar::computeExactError(Level& level, ConstVector<double> solution, Vector<double> error,
                                                       const ExactSolution& exact_solution)
{
    const PolarGrid& grid        = level.grid();
    const LevelCache& levelCache = level.levelCache();

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
