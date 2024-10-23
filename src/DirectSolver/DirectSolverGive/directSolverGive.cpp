#include "../../../include/DirectSolver/DirectSolverGive/directSolverGive.h"

// clang-format off
DirectSolverGive::DirectSolverGive(
    const PolarGrid& grid, const LevelCache& level_cache, 
    const DomainGeometry& domain_geometry, const DensityProfileCoefficients& density_profile_coefficients,
    bool DirBC_Interior, int num_omp_threads
) :
    DirectSolver(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads)
{
    solver_matrix_ = buildSolverMatrix();
    initializeMumpsSolver(mumps_solver_, solver_matrix_);
}

void DirectSolverGive::solveInPlace(Vector<double>& solution) 
{
    // Adjusts the right-hand side vector to account for symmetry corrections.
    // This transforms the system matrixA * solution = rhs into the equivalent system:
    // symmetric(matrixA) * solution = rhs - applySymmetryShift(rhs).
    applySymmetryShift(solution);
    // Solves the adjusted system symmetric(matrixA) * solution = rhs using the MUMPS solver.
    solveWithMumps(solution);
}


DirectSolverGive::~DirectSolverGive() 
{
    finalizeMumpsSolver(mumps_solver_);
}