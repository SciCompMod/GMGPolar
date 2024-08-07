#include "../../include/GMGPolar/gmgpolar.h"

#include <chrono>

void GMGPolar::solve() {
    auto start_solve = std::chrono::high_resolution_clock::now();

    int start_level_depth = 0;
    Level& level = levels_[start_level_depth];

    /* ---------------------------- */
    /* Initialize starting solution */
    assign(level.solution(), 0.0);

// Iteration: 0
// Exact Weighted-Euclidean Error: 3.2363237151e-05
// Exact Infinity Error: 9.9999999964e-05
// Residual Norm: 1.2822691453e-03
// Relative Residual Norm: 1.0000000000e+00

// Iteration: 0
// Exact Weighted-Euclidean Error: 3.2363237151e-05
// Exact Infinity Error: 9.9999999964e-05
// Residual Norm: 1.2822691453e-03
// Relative Residual Norm: 1.0000000000e+00

    {
        assign(level.solution(), 1.0);

    auto& grid = level.grid();
    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++)
    {
        int i_r;
        double r, theta, sin_theta, cos_theta;

        i_r = 0;
        r = grid.radius(i_r);
        theta = grid.theta(i_theta);
        sin_theta = sin(theta);
        cos_theta = cos(theta);
        level.solution()[grid.index(i_r, i_theta)] = exact_solution_->exact_solution(r, theta, sin_theta, cos_theta);

        i_r = grid.nr()-1;
        r = grid.radius(i_r);
        theta = grid.theta(i_theta);
        sin_theta = sin(theta);
        cos_theta = cos(theta);
        level.solution()[grid.index(i_r, i_theta)] = exact_solution_->exact_solution(r, theta, sin_theta, cos_theta);
    }
    }

    std::cout<<level.solution()<<std::endl;

    /* ---------------- */
    /* Discretize rhs_f */
    if(extrapolation_ > 0){
        Level& next_level = levels_[start_level_depth+1];
        injectToLowerLevel(start_level_depth, next_level.rhs(), level.rhs());
        discretize_rhs_f(next_level, next_level.rhs());
    }
    discretize_rhs_f(level, level.rhs());

    /* ------------ */
    /* Start Solver */

    int number_of_iterations_ = 0;

    double initial_residual_norm; 
    double current_residual_norm, current_relative_residual_norm;

    while(number_of_iterations_ < max_iterations_){
        std::cout<<"\nIteration: "<< number_of_iterations_ << std::endl;

        /* ---------------------------------------------- */
        /* Test solution against exact solution if given. */
        /* ---------------------------------------------- */
        if(exact_solution_ != nullptr) {
            auto start_check_exact_error = std::chrono::high_resolution_clock::now();

            std::pair<double, double> exact_error = compute_exact_error(level, level.solution(), level.residual());
            exact_errors_.push_back(exact_error);

            auto end_check_exact_error = std::chrono::high_resolution_clock::now();
            t_check_exact_error += std::chrono::duration<double>(end_check_exact_error - start_check_exact_error).count();

            std::cout << "Exact Weighted-Euclidean Error: " << exact_error.first << std::endl;
            std::cout << "Exact Infinity Error: " << exact_error.second << std::endl;
        }

        /* ---------------------------- */
        /* Compute convergence criteria */
        /* ---------------------------- */
        if(absolute_tolerance_.has_value() || relative_tolerance_.has_value()){
            auto start_check_convergence = std::chrono::high_resolution_clock::now();

            level.computeResidual(level.residual(), level.rhs(), level.solution());
            if(extrapolation_){
                Level& next_level = levels_[start_level_depth+1];
                injectToLowerLevel(start_level_depth, next_level.solution(), level.solution());
                next_level.computeResidual(next_level.residual(), next_level.rhs(), next_level.solution());
                extrapolated_residual(start_level_depth, level.residual(), next_level.residual());
            }

            switch (residual_norm_type_)
            {
                case ResidualNormType::EUCLIDEAN:
                    current_residual_norm = sqrt(l2_norm_squared(level.residual()));
                    break;
                case ResidualNormType::WEIGHTED_EUCLIDEAN:
                    current_residual_norm = sqrt(l2_norm_squared(level.residual())) / sqrt(level.grid().number_of_nodes());
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
            } else{
                current_relative_residual_norm = current_residual_norm / initial_residual_norm;
            }

            auto end_check_convergence = std::chrono::high_resolution_clock::now();
            t_check_convergence += std::chrono::duration<double>(end_check_convergence - start_check_convergence).count();

            std::cout<< "Residual Norm: " << current_residual_norm <<std::endl;
            std::cout<< "Relative Residual Norm: " << current_relative_residual_norm <<std::endl;

            if(converged(current_residual_norm, current_relative_residual_norm)) break;
        }

        /* ------------------------- */
        /* Start Multigrid Iteration */
        /* ------------------------- */
        auto start_solve_multigrid_iterations = std::chrono::high_resolution_clock::now();

        switch (multigrid_cycle_)
        {
            case MultigridCycleType::V_CYCLE:
                if(!extrapolation_){
                    multigrid_V_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                } else{
                    implicitly_extrapolated_multigrid_V_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                }
                break;
            case MultigridCycleType::W_CYCLE:
                if(!extrapolation_){
                    multigrid_W_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                } else{
                    implicitly_extrapolated_multigrid_W_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                }
                break;
            case MultigridCycleType::F_CYCLE:
                if(!extrapolation_){
                    multigrid_F_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                } else{
                    implicitly_extrapolated_multigrid_F_Cycle(start_level_depth, level.solution(), level.rhs(), level.residual());
                }
                break;
            default:
                throw std::invalid_argument("Unknown MultigridCycleType");
        }
        number_of_iterations_++;

        // std::cout<<"HELLO"<<std::endl;
        // std::cout<<level.solution()<<std::endl;

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
        mean_residual_reduction_factor_rho_ = std::pow(current_residual_norm / initial_residual_norm, 1.0 / number_of_iterations_);

        std::cout<<"\nTotal Iterations: "<<number_of_iterations_<<std::endl;
        std::cout<<"Mean Residual Reduction Factor Rho: "<< mean_residual_reduction_factor_rho_ <<std::endl;
    }

    auto end_solve = std::chrono::high_resolution_clock::now();
    t_solve_total += std::chrono::duration<double>(end_solve - start_solve).count();
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


std::pair<double, double> GMGPolar::compute_exact_error(Level& level, const Vector<double>& solution, Vector<double>& error){
    const PolarGrid& grid = level.grid();
    const LevelCache& levelCache = level.levelCache();
    const auto& sin_theta_cache = levelCache.sin_theta();
    const auto& cos_theta_cache = levelCache.cos_theta();

    assert(solution.size() == error.size());
    assert(solution.size() == grid.number_of_nodes());

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

    double weighted_euclidean_error = sqrt(l2_norm_squared(error)) / sqrt(grid.number_of_nodes());
    double infinity_error = infinity_norm(error);

    return std::make_pair(weighted_euclidean_error, infinity_error);
}


void GMGPolar::extrapolated_residual(const int current_level, Vector<double>& residual, const Vector<double>& residual_next_level){
    omp_set_num_threads(maxOpenMPThreads_);

    const PolarGrid& fineGrid = levels_[current_level].grid();
    const PolarGrid& coarseGrid = levels_[current_level+1].grid();

    assert(residual.size() == fineGrid.number_of_nodes());
    assert(residual_next_level.size() == coarseGrid.number_of_nodes());

    #pragma omp parallel num_threads(maxOpenMPThreads_)
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
