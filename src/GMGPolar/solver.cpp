#include "../../include/GMGPolar/gmgpolar.h"


#include <chrono>
#include <iostream>

// Assuming the other necessary includes and class definitions are already present

// void GMGPolar::multigrid_iteration(const int level_depth) {
//     assert(0 <= level_depth && level_depth < numberOflevels_-1);

//     Level& level = levels_[level_depth];
//     Level& next_level = levels_[level_depth+1];

//     auto start = std::chrono::high_resolution_clock::now();

//     /* ------------ */
//     /* Presmoothing */
//     auto presmooth_start = std::chrono::high_resolution_clock::now();
//     for (int i = 0; i < v1; i++){
//         level.smoothingInPlace(level.solution(), level.levelCache().rhs(), level.residual());
//     }
//     auto presmooth_end = std::chrono::high_resolution_clock::now();

//     /* ---------------------- */
//     /* Coarse grid correction */
//     /* ---------------------- */

//     /* Compute the residual */
//     auto residual_start = std::chrono::high_resolution_clock::now();
//     level.computeResidual(level.residual(), level.levelCache().rhs(), level.solution());
//     // if(level_depth==0) std::cout<<  sqrt(dot_product(level.residual(), level.residual()) / level.residual().size())<<std::endl;
//     auto residual_end = std::chrono::high_resolution_clock::now();

//     /* Restrict the residual */
//     auto restrict_start = std::chrono::high_resolution_clock::now();
//     restrictToLowerLevel(level_depth, next_level.residual(), level.residual());
//     auto restrict_end = std::chrono::high_resolution_clock::now();

//     /* Solve A * error = residual */
//     auto solve_start = std::chrono::high_resolution_clock::now();
//     if(level_depth+1 == numberOflevels_-1){
//         /* Using a direct solver */
//         next_level.coarseSolveInPlace(next_level.residual());
//     } else{
//         /* Recursively calling the multigrid cycle */
//         next_level.levelCache().rhs() = next_level.residual();
//         assign(next_level.solution(), 0.0);
//         multigrid_iteration(level_depth+1);
//     }
//     auto solve_end = std::chrono::high_resolution_clock::now();

//     /* Interpolate the correction */
//     auto prolongate_start = std::chrono::high_resolution_clock::now();
//     prolongateToNextLevel(level_depth+1, level.residual(), next_level.residual());
//     auto prolongate_end = std::chrono::high_resolution_clock::now();

//     /* Compute the corrected approximation: u = u + error */
//     auto add_start = std::chrono::high_resolution_clock::now();
//     add(level.solution(), level.residual());
//     auto add_end = std::chrono::high_resolution_clock::now();

//     /* ------------- */
//     /* Postsmoothing */
//     auto postsmooth_start = std::chrono::high_resolution_clock::now();
//     for (int i = 0; i < v2; i++){
//         level.smoothingInPlace(level.solution(), level.levelCache().rhs(), level.residual());
//     }
//     auto postsmooth_end = std::chrono::high_resolution_clock::now();

//     auto end = std::chrono::high_resolution_clock::now();

//     std::cout<<level_depth<<std::endl;

//     // Print durations
//     std::cout << "Presmoothing time: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(presmooth_end - presmooth_start).count() / 1000000
//               << " seconds" << std::endl;
//     std::cout << "Residual computation time: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(residual_end - residual_start).count() / 1000000
//               << " seconds" << std::endl;
//     std::cout << "Restriction time: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(restrict_end - restrict_start).count() / 1000000
//               << " seconds" << std::endl;
//     std::cout << "Solving time: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(solve_end - solve_start).count() / 1000000
//               << " seconds" << std::endl;
//     std::cout << "Prolongation time: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(prolongate_end - prolongate_start).count() / 1000000
//               << " seconds" << std::endl;
//     std::cout << "Addition time: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(add_end - add_start).count() / 1000000
//               << " seconds" << std::endl;
//     std::cout << "Postsmoothing time: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(postsmooth_end - postsmooth_start).count() / 1000000
//               << " seconds" << std::endl;
//     std::cout << "Total time: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000
//               << " seconds" << std::endl;
// }


void GMGPolar::multigrid_iteration(const int level_depth) {
    assert(0 <= level_depth && level_depth < numberOflevels_-1);

    Level& level = levels_[level_depth];
    Level& next_level = levels_[level_depth+1];

    /* ------------ */
    /* Presmoothing */
    for (int i = 0; i < v1; i++){
        level.smoothingInPlace(level.solution(), level.levelCache().rhs(), level.residual());
    }

    /* ---------------------- */
    /* Coarse grid correction */
    /* ---------------------- */

    /* Compute the residual */
    level.computeResidual(level.residual(), level.levelCache().rhs(), level.solution());

    /* Restrict the residual */
    restrictToLowerLevel(level_depth, next_level.residual(), level.residual());

    /* Solve A * error = residual */
    if(level_depth+1 == numberOflevels_-1){
        /* Using a direct solver */
        next_level.coarseSolveInPlace(next_level.residual());
    } else{
        /* Recursively calling the multigrid cycle */
        next_level.levelCache().rhs() = next_level.residual();
        assign(next_level.solution(), 0.0);
        multigrid_iteration(level_depth+1);
    }

    /* Interpolate the correction */
    prolongateToNextLevel(level_depth+1, level.residual(), next_level.residual());

    /* Compute the corrected approximation: u = u + error */
    add(level.solution(), level.residual());

    /* ------------- */
    /* Postsmoothing */
    for (int i = 0; i < v2; i++){
        level.smoothingInPlace(level.solution(), level.levelCache().rhs(), level.residual());
    }
}


void GMGPolar::solve() {

    int start_level_depth = 0;
    Level& level = levels_[start_level_depth];

    if(!extrapolation){

    
    }


    assign(level.solution(), 0.0);
    level.computeResidual(level.residual(), level.levelCache().rhs(), level.solution());
    std::cout<< sqrt(dot_product(level.residual(), level.residual())) / sqrt(level.grid().number_of_nodes()) <<std::endl;

    double factor = 1.0 / sqrt(level.grid().number_of_nodes());

    for (int i = 0; i < 100; i++)
    {
        multigrid_iteration(start_level_depth);

        level.computeResidual(level.residual(), level.levelCache().rhs(), level.solution());
        std::cout<< sqrt(dot_product(level.residual(), level.residual())) / sqrt(level.grid().number_of_nodes()) <<std::endl;

    }
    

    // write_to_vtk("SolutionMGC", level.grid(), level.solution(), level.levelCache());






    // int level_depth = 0;

    // const auto& grid = levels_[level_depth].grid();
    // const int n = grid.number_of_nodes();

    // Vector<double> error(n);

    // for (int i = 0; i < 1; i++)
    // {
    //     for (int index = 0; index < grid.number_of_nodes(); index++)
    //     {
    //         MultiIndex node = grid.multiindex(index);
    //         Point coords = grid.polar_coordinates(node);
    //         error[index] = std::abs( (*exact_solution_).exact_solution(coords[0], coords[1], sin(coords[1]), cos(coords[1])) - level.solution()[index]);

    //         error[index] = (*exact_solution_).exact_solution(coords[0], coords[1], sin(coords[1]), cos(coords[1]));
    //     }

    //     // #pragma omp parallel for
    //     // for (size_t i = 0; i < n; i++) error[i] = std::abs(solution[i] - y[i]);
    //     std::cout<<dot_product(error,error)<<std::endl;
    // }

    // write_to_vtk("errorMGC", grid, error, level.levelCache());














    // Level& level = levels_[level_depth];

    // Vector<double>& f_0 = level.levelCache().discretization_rhs_f();
    // Vector<double>& u_0 = level.solution();
    // Vector<double>& r_0 = level.residual();

    // assign(u_0, 0.0); 

    // multigrid_iteration(level_depth);





    // /* Presmoothing */
    // for (int i = 0; i < v1; i++){
    //     r_0 = f_0;
    //     level.smoothingInPlace(u_0, r_0);
    // }

    // /* Residual */
    // level.computeResidual(r_0, u_0);


    // Level& next_level = levels_[level_depth-1];
    // Vector<double>& temp_1 = level.levelCache().discretization_rhs_f();
    // Vector<double>& e_1 = level.solution();
    // Vector<double>& r_1 = level.residual();

    // /* Restrict Residual */
    // restrictToLowerLevel(level_depth, r_1, r_0);




    // /* Presmoothing */
    // for (int i = 0; i < v1; i++){
    //     r_1 = temp_1;
    //     next_level.smoothingInPlace(u_0, r_0);
    // }



    // int current_level_depth = 0;
    // Level& current_level = levels_[current_level_depth];

    // const auto& grid = current_level.grid();
    // const int n = grid.number_of_nodes();

    // Vector<double> solution(n);
    // Vector<double> temp(n);

    // assign(solution, 0.0);

    // /* Presmoothing */
    // for (int i = 0; i < v1; i++){
    //     current_level.smoothingInPlace(solution, temp);
    // }

    // current_level.computeResidual(temp, solution);


    // int next_level_depth = current_level_depth + 1;
    // Level& next_level = levels_[next_level_depth];

    // const auto& next_grid = next_level.grid();
    // const int next_n = next_grid.number_of_nodes();

    // Vector<double> solution2(next_n);
    // Vector<double> temp2(next_n);

    // restrictToLowerLevel(current_level_depth, temp2, temp);

    // coa



    // int current_level = 0;

    // const auto& grid = levels_[current_level].grid();
    // const int n = grid.number_of_nodes();

    // Vector<double> x(n);
    // Vector<double> result(levels_[current_level+1].grid().number_of_nodes());

    // assign(x, 4.0);

    // for (int i = 0; i < x.size(); i++)
    // {
    //     x[i] = 4.0 * i;
    // }

    

    // // interpolation_ -> applyRestrictionTake0(levels_[current_level], levels_[current_level+1], result, x);

    // restrictToLowerLevel(current_level, result, x);

    // interpolation_ -> applyProlongation0(levels_[current_level+1], levels_[current_level], x, result);

    // std::cout<< result<<std::endl;

    // std::cout<<x<<std::endl;

}

    // int current_level = 0;

    // const auto& grid = levels_[current_level].grid();
    // const int n = grid.number_of_nodes();

    // Vector<double> x(n);
    // Vector<double> y(n);

    // // Timing variables
    // auto start = std::chrono::high_resolution_clock::now();
    // auto end = start;

    // // Step 1: Assign x
    // start = std::chrono::high_resolution_clock::now();
    // assign(x, 0.0);
    // end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> time_assign = end - start;
    
    // // Step 2: Smoothing
    // start = std::chrono::high_resolution_clock::now();
    // levels_[current_level].smoothingInPlace(x, y);
    // end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> time_smoothing = end - start;

    // // Step 3: Compute Residual
    // start = std::chrono::high_resolution_clock::now();
    // levels_[current_level].computeResidual(y, x);
    // end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> time_compute_residual = end - start;

    // // Step 4: Coarse Solve
    // start = std::chrono::high_resolution_clock::now();
    // levels_[numberOflevels_-1].coarseSolveInPlace(y);
    // end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> time_coarse_solve = end - start;

    // // Step 5: Update x
    // start = std::chrono::high_resolution_clock::now();
    // #pragma omp parallel for
    // for (int i = 0; i < x.size(); i++) {
    //     x[i] += y[i];
    // }
    // end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> time_update_x = end - start;

    // // Print times
    // std::cout << "Time for assign: " << time_assign.count() << " seconds" << std::endl;
    // std::cout << "Time for smoothing: " << time_smoothing.count() << " seconds" << std::endl;
    // std::cout << "Time for compute residual: " << time_compute_residual.count() << " seconds" << std::endl;
    // std::cout << "Time for coarse solve: " << time_coarse_solve.count() << " seconds" << std::endl;
    // std::cout << "Time for update x: " << time_update_x.count() << " seconds" << std::endl;

    // // std::cout<<x[0]<<" "<< y[0]<<std::endl;

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

