#include "../../../include/DirectSolver/DirectSolverGiveCustomLU/directSolverGiveCustomLU.h"

DirectSolverGiveCustomLU::DirectSolverGiveCustomLU(const PolarGrid& grid, const LevelCache& level_cache,
                                                   const DomainGeometry& domain_geometry,
                                                   const DensityProfileCoefficients& density_profile_coefficients,
                                                   bool DirBC_Interior, int num_omp_threads)
    : DirectSolver(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads)
{
    solver_matrix_ = buildSolverMatrix();
    lu_solver_     = SparseLUSolver<double>(solver_matrix_);
}

void DirectSolverGiveCustomLU::solveInPlace(Vector<double>& solution)
{
    lu_solver_.solveInPlace(solution);
}

DirectSolverGiveCustomLU::~DirectSolverGiveCustomLU()
{
}