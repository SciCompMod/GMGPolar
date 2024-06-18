#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::solve() {






    int current_level = numberOflevels_ - 1;

    const auto& grid = levels_[current_level].grid();
    const int n = grid.number_of_nodes();

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
    // // levels_[current_level].coarseSolveInPlace(y);
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

    // std::cout<<x[0]<<" "<< y[0]<<std::endl;


    Vector<double> x(n);
    assign(x, 0.0);
    Vector<double> residual(n);

    levels_[current_level].computeResidual(residual,x);
    levels_[current_level].coarseSolveInPlace(residual);

    Vector<double> solution = residual;

    levels_[current_level].smoothingInPlace(solution, residual);



    Vector<double> y(n);
    assign(y, 0.0);
    
    y = solution;

    Vector<double> temp(n);

    Vector<double> error(n);

    for (int i = 0; i < 10; i++)
    {
        levels_[current_level].smoothingInPlace(y,temp);
        #pragma omp for
        for (size_t i = 0; i < n; i++) error[i] = std::abs(solution[i] - y[i]);
        std::cout<<dot_product(error,error)<<std::endl;
    }



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

}