#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::solve() {






    int current_level = 0;

    const auto& grid = levels_[current_level].grid();
    const int n = grid.number_of_nodes();

    Vector<double> x(n);
    assign(x, 0.0);
    Vector<double> residual(n);

    levels_[current_level].computeResidual(residual,x);
    levels_[current_level].coarseSolveInPlace(residual);

    Vector<double> solution = residual;


    Vector<double> y(n);
    assign(y, 0.0);
    
    // y = solution;

    Vector<double> temp(n);

    Vector<double> error(n);

    for (int i = 0; i < 300; i++)
    {
        levels_[current_level].smoothingInPlace(y,temp);
        #pragma omp for
        for (size_t i = 0; i < n; i++) error[i] = std::abs(solution[i] - y[i]);
        std::cout<<dot_product(error,error)<<std::endl;
    }


    // std::cout<< error<<std::endl;

    std::cout<<"Error Circle Section: "<<std::endl;

    for (int i = 0; i < grid.numberSmootherCircles(); i++)
    {
        for (int j = 0; j < grid.nr(); j++)
        {
            // std::cout<< error[grid.index(i,j)]<<std::endl;
        }
    }

    std::cout<<"Error Circle Section: "<<std::endl;

    for (int i = grid.numberSmootherCircles(); i < grid.nr(); i++)
    {
        for (int j = 0; j < grid.nr(); j++)
        {
            // std::cout<< error[grid.index(i,j)]<<std::endl;
        }
    }
    
    



    // levels_[current_level].smoothingInPlace(v2,v1);
    // levels_[current_level].smoothingInPlace(v2,v1);
    // levels_[current_level].smoothingInPlace(v2,v1);
    // levels_[current_level].smoothingInPlace(v2,v1);

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