#include "../../include/GMGPolar/gmgpolar.h"

#include <chrono>

void GMGPolar::solve() {
    auto start_solve = std::chrono::high_resolution_clock::now();

    /* ---------------- */
    /* Discretize rhs_f */
    /* ---------------- */

    int initial_rhs_f_levels = FMG_ ? number_of_levels_ : (extrapolation_ ? 2 : 1);
    // Loop through the levels, injecting and discretizing rhs
    for (int level_idx = 0; level_idx < initial_rhs_f_levels; ++level_idx) 
    {
        Level& current_level = levels_[level_idx];

        // Inject rhs if there is a next level
        if (level_idx + 1 < initial_rhs_f_levels) {
            Level& next_level = levels_[level_idx + 1];
            injection(level_idx, next_level.rhs(), current_level.rhs());
        }

        // Discretize the rhs for the current level
        discretize_rhs_f(current_level, current_level.rhs());
    }

    /* ---------------------------- */
    /* Initialize starting solution */
    /* ---------------------------- */

    if (!FMG_) {
        int start_level_depth = 0;
        Level& level = levels_[start_level_depth];
        assign(level.solution(), 0.0); // Assign zero initial guess if not using FMG
    }
    else {
        // Start from the coarsest level
        int FMG_start_level_depth = number_of_levels_ - 1;
        Level& FMG_level = levels_[FMG_start_level_depth];

        // Solve directly on the coarsest level
        FMG_level.solution() = FMG_level.rhs();
        FMG_level.directSolveInPlace(FMG_level.solution());  // Direct solve on coarsest grid

        // Prolongate the solution from the coarsest level up to the finest, while applying a V-Cycle on each level
        for (int current_level = FMG_start_level_depth - 1; current_level > 0; --current_level) {
            Level& FMG_level = levels_[current_level];         // The current level
            Level& next_FMG_level = levels_[current_level-1];  // The finer level

            prolongation(current_level, next_FMG_level.solution(), FMG_level.solution());
            multigrid_V_Cycle(current_level-1, next_FMG_level.solution(), next_FMG_level.rhs(), next_FMG_level.residual());
        }
    }

    /* ------------ */
    /* Start Solver */
    /* ------------ */

    int start_level_depth = 0;
    Level& level = levels_[start_level_depth];

    int number_of_iterations_ = 0;

    double initial_residual_norm; 
    double current_residual_norm, current_relative_residual_norm;

    while(number_of_iterations_ < max_iterations_){

        if(verbose_ > 0) {
            std::cout <<"\nIteration: "<< number_of_iterations_ << std::endl;
        }

        /* ---------------------------------------------- */
        /* Test solution against exact solution if given. */
        /* ---------------------------------------------- */
        if(exact_solution_ != nullptr) {
            auto start_check_exact_error = std::chrono::high_resolution_clock::now();

            std::pair<double, double> exact_error = computeExactError(level, level.solution(), level.residual());
            exact_errors_.push_back(exact_error);

            auto end_check_exact_error = std::chrono::high_resolution_clock::now();
            t_check_exact_error += std::chrono::duration<double>(end_check_exact_error - start_check_exact_error).count();

            if(verbose_ > 0) {
                std::cout << "Exact Weighted-Euclidean Error: " << exact_error.first << std::endl;
                std::cout << "Exact Infinity Error: " << exact_error.second << std::endl;
            }
        }

        /* ---------------------------- */
        /* Compute convergence criteria */
        /* ---------------------------- */
        if(absolute_tolerance_.has_value() || relative_tolerance_.has_value()){
            auto start_check_convergence = std::chrono::high_resolution_clock::now();

            level.computeResidual(level.residual(), level.rhs(), level.solution());
            if(extrapolation_){
                Level& next_level = levels_[start_level_depth+1];
                injection(start_level_depth, next_level.solution(), level.solution());
                next_level.computeResidual(next_level.residual(), next_level.rhs(), next_level.solution());
                extrapolatedResidual(start_level_depth, level.residual(), next_level.residual());
            }

            switch (residual_norm_type_)
            {
                case ResidualNormType::EUCLIDEAN:
                    current_residual_norm = sqrt(l2_norm_squared(level.residual()));
                    break;
                case ResidualNormType::WEIGHTED_EUCLIDEAN:
                    current_residual_norm = sqrt(l2_norm_squared(level.residual())) / sqrt(level.grid().numberOfNodes());
                    break;
                case ResidualNormType::INFINITY_NORM:
                    current_residual_norm = infinity_norm(level.residual());
                    break;
                default:
                    throw std::invalid_argument("Unknown ResidualNormType");
            }
            residual_norms_.push_back(current_residual_norm);
            
            if(number_of_iterations_ == 0) {
                initial_residual_norm = current_residual_norm;
                current_relative_residual_norm = 1.0;
                if(verbose_ > 0) {
                    std::cout << "Residual Norm: " << current_residual_norm << std::endl;
                }
            } 
            else{
                current_relative_residual_norm = current_residual_norm / initial_residual_norm;
                const double current_residual_reduction_factor = residual_norms_[number_of_iterations_] / residual_norms_[number_of_iterations_-1];

                if(verbose_ > 0) {
                    std::cout << "Residual Norm: " << current_residual_norm << std::endl;
                    std::cout << "Relative Residual Norm: " << current_relative_residual_norm << std::endl;
                    std::cout << "Residual Reduction Factor: " << current_residual_reduction_factor << std::endl;
                }

                bool switched_extrapolation = false;
                if(current_residual_reduction_factor > 0.4 && extrapolation_ == 3){
                    switched_extrapolation = true;
                    extrapolation_ = 1;
                    if(verbose_ > 0){
                        std::cout << "Switching from full grid smoothing to standard extrapolated smoothing!" << std::endl;
                    }
                }

                if(!switched_extrapolation && (current_residual_reduction_factor > 1.0 || equals(current_residual_reduction_factor, 1.0))){
                    if(verbose_ > 0){
                        std::cout << "Stopped Multigrid: Residuals are not improving further!" << std::endl;
                    }
                    break;
                }
            }

            auto end_check_convergence = std::chrono::high_resolution_clock::now();
            t_check_convergence += std::chrono::duration<double>(end_check_convergence - start_check_convergence).count();

            if(converged(current_residual_norm, current_relative_residual_norm)) break;
        }

        /* ------------------------- */
        /* Start Multigrid Iteration */
        /* ------------------------- */
        auto start_solve_multigrid_iterations = std::chrono::high_resolution_clock::now();

        switch (multigrid_cycle_)
        {
            case MultigridCycleType::V_CYCLE:
                if(extrapolation_ == 0){
                    multigrid_V_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                } else{
                    implicitlyExtrapolatedMultigrid_V_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                }
                break;
            case MultigridCycleType::W_CYCLE:
                if(extrapolation_ == 0){
                    multigrid_W_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                } else{
                    implicitlyExtrapolatedMultigrid_W_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                }
                break;
            case MultigridCycleType::F_CYCLE:
                if(extrapolation_ == 0){
                    multigrid_F_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                } else{
                    implicitlyExtrapolatedMultigrid_F_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                }
                break;
            default:
                throw std::invalid_argument("Unknown MultigridCycleType");
        }
        number_of_iterations_++;

        auto end_solve_multigrid_iterations = std::chrono::high_resolution_clock::now();
        t_solve_multigrid_iterations += std::chrono::duration<double>(end_solve_multigrid_iterations - start_solve_multigrid_iterations).count();
    }

    if(number_of_iterations_ > 0){
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
        mean_residual_reduction_factor_ = std::pow(current_residual_norm / initial_residual_norm, 1.0 / number_of_iterations_);

        if(verbose_ > 0) {
            std::cout <<"\nTotal Iterations: "<<number_of_iterations_<< std::endl;
            std::cout <<"Mean Residual Reduction Factor Rho: "<< mean_residual_reduction_factor_ << std::endl;
        }
    }

    auto end_solve = std::chrono::high_resolution_clock::now();
    t_solve_total += std::chrono::duration<double>(end_solve - start_solve).count();

    if(verbose_ > 50){
        computeExactError(level, level.solution(), level.residual());
        writeToVTK("solution", level, level.solution());
        writeToVTK("error", level, level.residual());
    }
}

 
bool GMGPolar::converged(const double& residual_norm, const double& relative_residual_norm){
    if(relative_tolerance_.has_value()){
        if(!(relative_residual_norm > relative_tolerance_.value())){
            return true;
        }
    }
    if(absolute_tolerance_.has_value()){
        if(!(residual_norm > absolute_tolerance_.value())){
            return true;
        }
    }
    return false;
}


std::pair<double, double> GMGPolar::computeExactError(Level& level, const Vector<double>& solution, Vector<double>& error){
    assert(exact_solution_ != nullptr);

    omp_set_num_threads(threads_per_level_[level.level()]);

    const PolarGrid& grid = level.grid();
    const LevelCache& levelCache = level.levelCache();
    const auto& sin_theta_cache = levelCache.sin_theta();
    const auto& cos_theta_cache = levelCache.cos_theta();

    assert(solution.size() == error.size());
    assert(solution.size() == grid.numberOfNodes());

    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++){
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                double theta = grid.theta(i_theta);
                double sin_theta = sin_theta_cache[i_theta];
                double cos_theta = cos_theta_cache[i_theta];
                error[grid.index(i_r, i_theta)] = exact_solution_->exact_solution(r, theta, sin_theta, cos_theta) - solution[grid.index(i_r, i_theta)];
            }
        }
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
            double theta = grid.theta(i_theta);
            double sin_theta = sin_theta_cache[i_theta];
            double cos_theta = cos_theta_cache[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                double r = grid.radius(i_r);
                error[grid.index(i_r, i_theta)] = exact_solution_->exact_solution(r, theta, sin_theta, cos_theta) - solution[grid.index(i_r, i_theta)];
            }
        }
    }

    double weighted_euclidean_error = sqrt(l2_norm_squared(error)) / sqrt(grid.numberOfNodes());
    double infinity_error = infinity_norm(error);

    return std::make_pair(weighted_euclidean_error, infinity_error);
}


void GMGPolar::extrapolatedResidual(const int current_level, Vector<double>& residual, const Vector<double>& residual_next_level){
    omp_set_num_threads(threads_per_level_[current_level]);

    const PolarGrid& fineGrid = levels_[current_level].grid();
    const PolarGrid& coarseGrid = levels_[current_level+1].grid();

    assert(residual.size() == fineGrid.numberOfNodes());
    assert(residual_next_level.size() == coarseGrid.numberOfNodes());

    #pragma omp parallel
    {
        /* Circluar Indexing Section */
        /* For loop matches circular access pattern */
        #pragma omp for nowait
        for (int i_r = 0; i_r < fineGrid.numberSmootherCircles(); i_r++){
            int i_r_coarse = i_r >> 1;
            for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
                int i_theta_coarse = i_theta >> 1;

                if(i_r & 1 || i_theta & 1){
                    residual[fineGrid.index(i_r,i_theta)] *= 4.0 / 3.0;
                }
                else{
                    int fine_idx = fineGrid.index(i_r, i_theta);
                    int coarse_idx = coarseGrid.index(i_r_coarse, i_theta_coarse);
                    residual[fine_idx] = (4.0 * residual[fine_idx] - residual_next_level[coarse_idx]) / 3.0;
                }   
            }
        }

        /* Radial Indexing Section */
        /* For loop matches radial access pattern */
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
            int i_theta_coarse = i_theta >> 1;
            for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                int i_r_coarse = i_r >> 1;

                if(i_r & 1 || i_theta & 1){
                    residual[fineGrid.index(i_r,i_theta)] *= 4.0 / 3.0;
                }
                else{
                    int fine_idx = fineGrid.index(i_r, i_theta);
                    int coarse_idx = coarseGrid.index(i_r_coarse, i_theta_coarse);
                    residual[fine_idx] = (4.0 * residual[fine_idx] - residual_next_level[coarse_idx]) / 3.0;
                }   
            }
        }
    }
}
