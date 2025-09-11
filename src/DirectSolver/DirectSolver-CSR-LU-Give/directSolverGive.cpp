#include "../../../include/DirectSolver/DirectSolver-CSR-LU-Give/directSolverGiveCustomLU.h"

DirectSolverGiveCustomLU::DirectSolverGiveCustomLU(const PolarGrid& grid, const LevelCache& level_cache,
                                                   const DomainGeometry& domain_geometry,
                                                   const DensityProfileCoefficients& density_profile_coefficients,
                                                   bool DirBC_Interior, int num_omp_threads)
    : DirectSolver(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads)
{
    solver_matrix_ = buildSolverMatrix();
    lu_solver_     = SparseLUSolver<double>(solver_matrix_);
}

void DirectSolverGiveCustomLU::solveInPlace(Vector<double> solution)
{
            std::cout<<" ENTER in SOLVEINPLACE  DirectSolverGiveCustomLU "<<std::endl;
 std::cout<<" solution(0) before  "<< solution(0)<<std::endl;
 std::cout<<" solution(1) before "<< solution(1)<<std::endl;


    lu_solver_.solveInPlace(solution);
    std::cout<<" solution(0) AFTER  "<< solution(0)<<std::endl;
 std::cout<<" solution(1) AFTER "<< solution(1)<<std::endl;
}

DirectSolverGiveCustomLU::~DirectSolverGiveCustomLU()
{
}
