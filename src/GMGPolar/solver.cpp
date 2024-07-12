#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::multigrid_V_Cycle(const int level_depth) {
    assert(0 <= level_depth && level_depth < numberOflevels_-1);

    auto start_MGC = std::chrono::high_resolution_clock::now();

    Level& level = levels_[level_depth];
    Level& next_level = levels_[level_depth+1];

    auto start_MGC_preSmoothing = std::chrono::high_resolution_clock::now();

    /* ------------ */
    /* Presmoothing */
    for (int i = 0; i < preSmoothingSteps_; i++){
        level.smoothingInPlace(level.solution(), level.rhs(), level.residual());
    }

    auto end_MGC_preSmoothing = std::chrono::high_resolution_clock::now();
    t_avg_MGC_preSmoothing += std::chrono::duration<double>(end_MGC_preSmoothing - start_MGC_preSmoothing).count();

    /* ---------------------- */
    /* Coarse grid correction */
    /* ---------------------- */

    auto start_MGC_residual = std::chrono::high_resolution_clock::now();

    /* Compute the residual */
    level.computeResidual(level.residual(), level.rhs(), level.solution());

    auto end_MGC_residual = std::chrono::high_resolution_clock::now();
    t_avg_MGC_residual += std::chrono::duration<double>(end_MGC_residual - start_MGC_residual).count();

    /* Restrict the residual */
    restrictToLowerLevel(level_depth, next_level.residual(), level.residual());

    /* Solve A * error = residual */
    if(level_depth+1 == numberOflevels_-1){
        auto start_MGC_directSolver = std::chrono::high_resolution_clock::now();

        /* Using a direct solver */
        next_level.directSolveInPlace(next_level.residual());

        auto end_MGC_directSolver = std::chrono::high_resolution_clock::now();
        t_avg_MGC_directSolver += std::chrono::duration<double>(end_MGC_directSolver - start_MGC_directSolver).count();
    } else{
        /* By recursively calling the multigrid cycle */
        next_level.rhs() = next_level.residual();
        assign(next_level.solution(), 0.0);
        multigrid_V_Cycle(level_depth+1);
    }

    /* Interpolate the correction */
    prolongateToUpperLevel(level_depth+1, level.residual(), next_level.residual());

    /* Compute the corrected approximation: u = u + error */
    add(level.solution(), level.residual());

    auto start_MGC_postSmoothing = std::chrono::high_resolution_clock::now();

    /* ------------- */
    /* Postsmoothing */
    for (int i = 0; i < postSmoothingSteps_; i++){
        level.smoothingInPlace(level.solution(), level.rhs(), level.residual());
    }

    auto end_MGC_postSmoothing = std::chrono::high_resolution_clock::now();
    t_avg_MGC_postSmoothing += std::chrono::duration<double>(end_MGC_postSmoothing - start_MGC_postSmoothing).count();

    auto end_MGC = std::chrono::high_resolution_clock::now();
    t_avg_MGC_total += std::chrono::duration<double>(end_MGC - start_MGC).count();
}


void GMGPolar::solve() {
    auto start_solve = std::chrono::high_resolution_clock::now();

    int start_level_depth = 0;
    Level& level = levels_[start_level_depth];

    /* ---------------------------- */
    /* Initialize starting solution */
    assign(level.solution(), 0.0);

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
                multigrid_V_Cycle(start_level_depth);
                break;
            case MultigridCycleType::W_CYCLE:
                break;
            case MultigridCycleType::F_CYCLE:
                break;
            default:
                throw std::invalid_argument("Unknown MultigridCycleType");
        }
        number_of_iterations_++;

        auto end_solve_multigrid_iterations = std::chrono::high_resolution_clock::now();
        t_solve_multigrid_iterations += std::chrono::duration<double>(end_solve_multigrid_iterations - start_solve_multigrid_iterations).count();
    }

    /* --------------------------------------------- */
    /* Compute the average Multigrid Iteration times */
    /* --------------------------------------------- */
    t_avg_MGC_total /= number_of_iterations_;
    t_avg_MGC_preSmoothing /= number_of_iterations_;
    t_avg_MGC_postSmoothing /= number_of_iterations_;
    t_avg_MGC_residual /= number_of_iterations_;
    t_avg_MGC_directSolver /= number_of_iterations_;

    /* -------------------------------- */
    /* Compute the reduction factor rho */
    /* -------------------------------- */
    mean_residual_reduction_factor_rho_ = std::pow(current_residual_norm / initial_residual_norm, 1.0 / number_of_iterations_);

    std::cout<< "\nMean Residual Reduction Factor Rho: "<< mean_residual_reduction_factor_rho_ <<std::endl;

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


    // int current_level = 0;

    // const auto& grid = levels_[current_level].grid();
    // const int n = grid.number_of_nodes();

    // Vector<double> x(n);
    // assign(x, 0.0);
    // Vector<double> residual(n);

    // levels_[current_level].computeResidual(residual,x);
    // levels_[current_level].coarseSolveInPlace(residual);

    // Vector<double> solution = residual;

    // // levels_[current_level].smoothingInPlace(solution, residual);


    // Vector<double> y(n);
    // assign(y, 0.0);
    
    // y = solution;

    // y[10] = 900;

    // Vector<double> temp(n);

    // Vector<double> error(n);

    // for (int i = 0; i < 100; i++)
    // {
    //     for (int index = 0; index < grid.number_of_nodes(); index++)
    //     {
    //         MultiIndex node = grid.multiindex(index);
    //         Point coords = grid.polar_coordinates(node);
    //         //error[index] = std::abs(solution[index] - y[index]);
    //         error[index] = std::abs( (*exact_solution_).exact_solution(coords[0], coords[1], sin(coords[1]), cos(coords[1])) - y[index]);
    //     }
    //     temp = levels_[current_level].levelCache().discretization_rhs_f();
    //     levels_[current_level].smoothingInPlace(y,temp);   

    //     // #pragma omp parallel for
    //     // for (size_t i = 0; i < n; i++) error[i] = std::abs(solution[i] - y[i]);
    //     std::cout<<dot_product(error,error)<<std::endl;
    // }






    // for (int i_r = 0; i_r < grid.nr(); i_r++)
    // {
    //     for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++)
    //     {
    //         v2[grid.index(i_r,i_theta)] = std::abs(v2[grid.index(i_r,i_theta)] - (*exact_solution_).exact_solution( grid.radius(i_r), 
    //             grid.theta(i_theta), sin(grid.theta(i_theta)), cos(grid.theta(i_theta))  ) );
    //     }
    // }

    // std::cout<<v2<<std::endl;
    
    
    // const std::filesystem::path file_path = "ErrorDirBC";
    // const LevelCache& level_data = levels_[current_level].levelData();

    // write_to_vtk(file_path, grid, v2, level_data);


























    // // std::cout<<"End"<<std::endl;
    // // std::cout<<v2<<std::endl;

    // int i_r = grid.nr() / 3;
    // int i_theta = grid.ntheta() / 3;

    // std::cout<<v2[grid.index(i_r, i_theta)]<<std::endl;

    // std::cout<<(*exact_solution_).exact_solution(grid.radius(i_r), 
    //     grid.theta(i_theta), sin(grid.theta(i_theta)), cos(grid.theta(i_theta)))<<std::endl;;

