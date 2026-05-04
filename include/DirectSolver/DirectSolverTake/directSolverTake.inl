#pragma once

template <class LevelCacheType>
DirectSolverTake<LevelCacheType>::DirectSolverTake(const PolarGrid& grid, const LevelCacheType& level_cache,
                                                   bool DirBC_Interior, int num_omp_threads)
    : DirectSolver<LevelCacheType>(grid, level_cache, DirBC_Interior, num_omp_threads)
#ifdef GMGPOLAR_USE_MUMPS
    , system_solver_(buildSolverMatrix())
#else
    , system_matrix_(buildSolverMatrix())
    , system_solver_(system_matrix_)
#endif
{
}

template <class LevelCacheType>
void DirectSolverTake<LevelCacheType>::solveInPlace(Vector<double> solution)
{
    // Adjusts the right-hand side vector to account for symmetry corrections.
    // This transforms the system matrixA * solution = rhs into the equivalent system:
    // symmetric_DBc(matrixA) * solution = rhs - applySymmetryShift(rhs).
    // The correction modifies the rhs to account for the influence of the Dirichlet boundary conditions,
    // ensuring that the solution at the boundary is correctly adjusted and maintains the required symmetry.
    applySymmetryShift(solution);
    // Solves the adjusted system symmetric(matrixA) * solution = rhs using the MUMPS solver.
    system_solver_.solveInPlace(solution);
}
