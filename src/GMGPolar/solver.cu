#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/LinearAlgebra/Vector/vector_operations.h"
#include "../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

#include <chrono>

void GMGPolar::solve()
{
    auto start_solve = std::chrono::high_resolution_clock::now();

    /* ---------------------------- */
    /* Initialize starting solution */
    /* ---------------------------- */

    auto start_initial_approximation = std::chrono::high_resolution_clock::now();

    if (!FMG_) {
        int start_level_depth = 0;
        Level& level          = levels_[start_level_depth];
        if(level.processingType() == ProcessingType::GPU){
            assign(level.GPU_solution(), 0.0);
        } else{
            assign(level.solution(), 0.0);
        }
    }
    else {
        // Start from the coarsest level
        int FMG_start_level_depth = number_of_levels_ - 1;
        Level& FMG_level          = levels_[FMG_start_level_depth];
        assert(FMG_level.processingType() != ProcessingType::GPU);

        // Solve directly on the coarsest level
        FMG_level.solution() = FMG_level.rhs();
        FMG_level.directSolveInPlace(FMG_level.solution()); // Direct solve on coarsest grid

        // Prolongate the solution from the coarsest level up to the finest, while applying Multigrid Cycles on each level
        for (int current_level = FMG_start_level_depth - 1; current_level > 0; --current_level) {
            Level& FMG_level      = levels_[current_level]; // The current level
            Level& next_FMG_level = levels_[current_level - 1]; // The finer level

            if(FMG_level.processingType() == ProcessingType::CPU_HYBRID){
                copyHostToDevice(FMG_level.solution(), FMG_level.GPU_solution());
            }
            if(FMG_level.processingType() != ProcessingType::CPU){
                FMGInterpolation(current_level, next_FMG_level.GPU_solution(), FMG_level.GPU_solution());
            } else{
                FMGInterpolation(current_level, next_FMG_level.solution(), FMG_level.solution());
            }

            // Apply some FMG iterations
            for (int i = 0; i < FMG_iterations_; i++) {
                if (current_level - 1 == 0 && (extrapolation_ != ExtrapolationType::NONE)) {
                    switch (FMG_cycle_) {
                    case MultigridCycleType::V_CYCLE:
                        if(next_FMG_level.processingType() == ProcessingType::GPU){
                            implicitlyExtrapolatedMultigrid_V_Cycle(
                                current_level - 1, next_FMG_level.GPU_solution(), next_FMG_level.GPU_rhs(), next_FMG_level.GPU_residual());
                        } else {
                            implicitlyExtrapolatedMultigrid_V_Cycle(
                                current_level - 1, next_FMG_level.solution(), next_FMG_level.rhs(), next_FMG_level.residual());
                        }
                        break;
                    case MultigridCycleType::W_CYCLE:
                         if(next_FMG_level.processingType() == ProcessingType::GPU){
                            implicitlyExtrapolatedMultigrid_W_Cycle(
                                current_level - 1, next_FMG_level.GPU_solution(), next_FMG_level.GPU_rhs(), next_FMG_level.GPU_residual());
                        } else {
                            implicitlyExtrapolatedMultigrid_F_Cycle(
                                current_level - 1, next_FMG_level.solution(), next_FMG_level.rhs(), next_FMG_level.residual());
                        }
                        break;
                    case MultigridCycleType::F_CYCLE:
                         if(next_FMG_level.processingType() == ProcessingType::GPU){
                            implicitlyExtrapolatedMultigrid_F_Cycle(
                                current_level - 1, next_FMG_level.GPU_solution(), next_FMG_level.GPU_rhs(), next_FMG_level.GPU_residual());
                        } else {
                            implicitlyExtrapolatedMultigrid_F_Cycle(
                                current_level - 1, next_FMG_level.solution(), next_FMG_level.rhs(), next_FMG_level.residual());
                        }
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
                         if(next_FMG_level.processingType() == ProcessingType::GPU){
                            multigrid_V_Cycle(
                                current_level - 1, next_FMG_level.GPU_solution(), next_FMG_level.GPU_rhs(), next_FMG_level.GPU_residual());
                        } else {
                            multigrid_V_Cycle(
                                current_level - 1, next_FMG_level.solution(), next_FMG_level.rhs(), next_FMG_level.residual());
                        }
                        break;
                    case MultigridCycleType::W_CYCLE:
                         if(next_FMG_level.processingType() == ProcessingType::GPU){
                            multigrid_W_Cycle(
                                current_level - 1, next_FMG_level.GPU_solution(), next_FMG_level.GPU_rhs(), next_FMG_level.GPU_residual());
                        } else {
                            multigrid_W_Cycle(
                                current_level - 1, next_FMG_level.solution(), next_FMG_level.rhs(), next_FMG_level.residual());
                        }
                        break;
                    case MultigridCycleType::F_CYCLE:
                         if(next_FMG_level.processingType() == ProcessingType::GPU){
                            multigrid_F_Cycle(
                                current_level - 1, next_FMG_level.GPU_solution(), next_FMG_level.GPU_rhs(), next_FMG_level.GPU_residual());
                        } else {
                            multigrid_F_Cycle(
                                current_level - 1, next_FMG_level.solution(), next_FMG_level.rhs(), next_FMG_level.residual());
                        }
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

    /* These times are included in the initial approximation and don't count towards the multigrid cyclces. */
    t_avg_MGC_total         = 0.0;
    t_avg_MGC_preSmoothing  = 0.0;
    t_avg_MGC_postSmoothing = 0.0;
    t_avg_MGC_residual      = 0.0;
    t_avg_MGC_directSolver  = 0.0;

    auto end_initial_approximation = std::chrono::high_resolution_clock::now();
    t_solve_initial_approximation +=
        std::chrono::duration<double>(end_initial_approximation - start_initial_approximation).count();

    /* ------------ */
    /* Start Solver */
    /* ------------ */

    int start_level_depth = 0;
    Level& level          = levels_[start_level_depth];

    number_of_iterations_ = 0;

    double initial_residual_norm;
    double current_residual_norm, current_relative_residual_norm;

    while (number_of_iterations_ < max_iterations_) {

        if (verbose_ > 0) {
            std::cout << "\nIteration: " << number_of_iterations_ << std::endl;
        }

        /* ---------------------------------------------- */
        /* Test solution against exact solution if given. */
        /* ---------------------------------------------- */
        if (exact_solution_ != nullptr) {
            auto start_check_exact_error = std::chrono::high_resolution_clock::now();

            std::pair<double, double> exact_error;
            if(level.processingType() == ProcessingType::GPU){
                exact_error = computeExactError(level, level.GPU_solution(), level.GPU_residual());
            } else {
                exact_error = computeExactError(level, level.solution(), level.residual());
            }
            exact_errors_.push_back(exact_error);

            auto end_check_exact_error = std::chrono::high_resolution_clock::now();
            t_check_exact_error +=
                std::chrono::duration<double>(end_check_exact_error - start_check_exact_error).count();

            if (verbose_ > 0) {
                std::cout << "Exact Weighted-Euclidean Error: " << exact_error.first << std::endl;
                std::cout << "Exact Infinity Error: " << exact_error.second << std::endl;
            }
        }

        /* ---------------------------- */
        /* Compute convergence criteria */
        /* ---------------------------- */
        if (absolute_tolerance_.has_value() || relative_tolerance_.has_value()) {
            auto start_check_convergence = std::chrono::high_resolution_clock::now();

            if(level.processingType() == ProcessingType::GPU){
                level.computeResidual(level.GPU_residual(), level.GPU_rhs(), level.GPU_solution());
            } else {
                level.computeResidual(level.residual(), level.rhs(), level.solution());
            }

            if (extrapolation_ != ExtrapolationType::NONE) {
                Level& next_level = levels_[start_level_depth + 1];

                if(level.processingType() == ProcessingType::GPU){
                    assert(next_level.processingType() != ProcessingType::CPU);
                    injection(start_level_depth, next_level.GPU_solution(), level.GPU_solution());
                } else{
                    assert(next_level.processingType() == ProcessingType::CPU);
                    injection(start_level_depth, next_level.solution(), level.solution());
                }

                if(next_level.processingType() == ProcessingType::GPU){
                    next_level.computeResidual(next_level.GPU_residual(), next_level.GPU_rhs(), next_level.GPU_solution());
                } else{
                    next_level.computeResidual(next_level.residual(), next_level.rhs(), next_level.solution());
                }

                if(next_level.processingType() == ProcessingType::CPU_HYBRID){
                    copyHostToDevice(next_level.residual(), next_level.GPU_residual());
                }

                if(level.processingType() == ProcessingType::GPU){
                    extrapolatedResidual(start_level_depth, level.GPU_residual(), next_level.GPU_residual());
                } else{
                    extrapolatedResidual(start_level_depth, level.residual(), next_level.residual());
                }
            }

            switch (residual_norm_type_) {
            case ResidualNormType::EUCLIDEAN:
                if(level.processingType() == ProcessingType::GPU){
                    current_residual_norm = l2_norm(level.GPU_residual());
                } else{
                    current_residual_norm = l2_norm(level.residual());
                }
                break;
            case ResidualNormType::WEIGHTED_EUCLIDEAN:
                if(level.processingType() == ProcessingType::GPU){
                    current_residual_norm = l2_norm(level.GPU_residual()) / sqrt(level.grid().numberOfNodes());
                } else{
                    current_residual_norm = l2_norm(level.GPU_residual()) / sqrt(level.grid().numberOfNodes());
                }
                break;
            case ResidualNormType::INFINITY_NORM:
                if(level.processingType() == ProcessingType::GPU){
                    current_residual_norm = infinity_norm(level.GPU_residual());
                } else{
                    current_residual_norm = infinity_norm(level.residual());
                }
                break;
            default:
                throw std::invalid_argument("Unknown ResidualNormType");
            }
            residual_norms_.push_back(current_residual_norm);

            if (number_of_iterations_ == 0) {
                initial_residual_norm          = current_residual_norm;
                current_relative_residual_norm = 1.0;
                if (verbose_ > 0) {
                    std::cout << "Residual Norm: " << current_residual_norm << std::endl;
                }
            }
            else {
                current_relative_residual_norm = current_residual_norm / initial_residual_norm;
                const double current_residual_reduction_factor =
                    residual_norms_[number_of_iterations_] / residual_norms_[number_of_iterations_ - 1];

                if (verbose_ > 0) {
                    std::cout << "Residual Norm: " << current_residual_norm << std::endl;
                    std::cout << "Relative Residual Norm: " << current_relative_residual_norm << std::endl;
                    std::cout << "Residual Reduction Factor: " << current_residual_reduction_factor << std::endl;
                }

                const double convergence_factor = 0.7;
                if (current_residual_reduction_factor > convergence_factor &&
                    extrapolation_ == ExtrapolationType::COMBINED && full_grid_smoothing_) {
                    full_grid_smoothing_ = false;
                    std::cout << "Switching from full grid smoothing to standard extrapolated smoothing." << std::endl;
                }
            }

            auto end_check_convergence = std::chrono::high_resolution_clock::now();
            t_check_convergence +=
                std::chrono::duration<double>(end_check_convergence - start_check_convergence).count();

            if (converged(current_residual_norm, current_relative_residual_norm))
                break;
        }

        /* ------------------------- */
        /* Start Multigrid Iteration */
        /* ------------------------- */
        auto start_solve_multigrid_iterations = std::chrono::high_resolution_clock::now();

        switch (multigrid_cycle_) {
        case MultigridCycleType::V_CYCLE:
            if (extrapolation_ == ExtrapolationType::NONE) {
                if(level.processingType() == ProcessingType::GPU){
                    multigrid_V_Cycle(
                        start_level_depth, level.GPU_solution(), level.GPU_rhs(), level.GPU_residual());
                } else {
                    multigrid_V_Cycle(
                        start_level_depth, level.solution(), level.rhs(), level.residual());
                }
            }
            else {
                if(level.processingType() == ProcessingType::GPU){
                    implicitlyExtrapolatedMultigrid_V_Cycle(
                        start_level_depth, level.GPU_solution(), level.GPU_rhs(), level.GPU_residual());
                } else {
                    implicitlyExtrapolatedMultigrid_V_Cycle(
                        start_level_depth, level.solution(), level.rhs(), level.residual());
                }
            }
            break;
        case MultigridCycleType::W_CYCLE:
            if (extrapolation_ == ExtrapolationType::NONE) {
                if(level.processingType() == ProcessingType::GPU){
                    multigrid_W_Cycle(
                        start_level_depth, level.GPU_solution(), level.GPU_rhs(), level.GPU_residual());
                } else {
                    multigrid_W_Cycle(
                        start_level_depth, level.solution(), level.rhs(), level.residual());
                }
            }
            else {
                if(level.processingType() == ProcessingType::GPU){
                    implicitlyExtrapolatedMultigrid_W_Cycle(
                        start_level_depth, level.GPU_solution(), level.GPU_rhs(), level.GPU_residual());
                } else {
                    implicitlyExtrapolatedMultigrid_W_Cycle(
                        start_level_depth, level.solution(), level.rhs(), level.residual());
                }
            }
            break;
        case MultigridCycleType::F_CYCLE:
            if (extrapolation_ == ExtrapolationType::NONE) {
                if(level.processingType() == ProcessingType::GPU){
                    multigrid_F_Cycle(
                        start_level_depth, level.GPU_solution(), level.GPU_rhs(), level.GPU_residual());
                } else {
                    multigrid_F_Cycle(
                        start_level_depth, level.solution(), level.rhs(), level.residual());
                }
            }
            else {
                if(level.processingType() == ProcessingType::GPU){
                    implicitlyExtrapolatedMultigrid_F_Cycle(
                        start_level_depth, level.GPU_solution(), level.GPU_rhs(), level.GPU_residual());
                } else {
                    implicitlyExtrapolatedMultigrid_F_Cycle(
                        start_level_depth, level.solution(), level.rhs(), level.residual());
                }
            }
            break;
        default:
            throw std::invalid_argument("Unknown MultigridCycleType");
        }
        number_of_iterations_++;

        auto end_solve_multigrid_iterations = std::chrono::high_resolution_clock::now();
        t_solve_multigrid_iterations +=
            std::chrono::duration<double>(end_solve_multigrid_iterations - start_solve_multigrid_iterations).count();
    }

    if (number_of_iterations_ > 0) {
        /* --------------------------------------------- */
        /* Compute the average Multigrid Iteration times */
        /* --------------------------------------------- */
        t_avg_MGC_total = t_solve_multigrid_iterations / number_of_iterations_;
        t_avg_MGC_preSmoothing /= number_of_iterations_;
        t_avg_MGC_postSmoothing /= number_of_iterations_;
        t_avg_MGC_residual /= number_of_iterations_;
        t_avg_MGC_directSolver /= number_of_iterations_;

        /* -------------------------------- */
        /* Compute the reduction factor rho */
        /* -------------------------------- */
        mean_residual_reduction_factor_ =
            std::pow(current_residual_norm / initial_residual_norm, 1.0 / number_of_iterations_);

        if (verbose_ > 0) {
            std::cout << "\nTotal Iterations: " << number_of_iterations_ << std::endl;
            std::cout << "Mean Residual Reduction Factor Rho: " << mean_residual_reduction_factor_ << std::endl;
        }
    }

    auto end_solve = std::chrono::high_resolution_clock::now();
    t_solve_total += std::chrono::duration<double>(end_solve - start_solve).count();
    t_solve_total -= t_check_exact_error;
}

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
