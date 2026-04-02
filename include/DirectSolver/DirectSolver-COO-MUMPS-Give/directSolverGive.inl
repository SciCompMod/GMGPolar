#pragma once

#ifdef GMGPOLAR_USE_MUMPS

template <class LevelCacheType>
DirectSolver_COO_MUMPS_Give<LevelCacheType>::DirectSolver_COO_MUMPS_Give(const PolarGrid& grid,
                                                                         const LevelCacheType& level_cache,
                                                                         bool DirBC_Interior, int num_omp_threads)
    : DirectSolver<LevelCacheType>(grid, level_cache, DirBC_Interior, num_omp_threads)
    , mumps_solver_(buildSolverMatrix())
{
}

template <class LevelCacheType>
void DirectSolver_COO_MUMPS_Give<LevelCacheType>::solveInPlace(Vector<double> solution)
{
    // Adjusts the right-hand side vector to account for symmetry corrections.
    // This transforms the system matrixA * solution = rhs into the equivalent system:
    // symmetric_DBc(matrixA) * solution = rhs - applySymmetryShift(rhs).
    // The correction modifies the rhs to account for the influence of the Dirichlet boundary conditions,
    // ensuring that the solution at the boundary is correctly adjusted and maintains the required symmetry.
    applySymmetryShift(solution);
    // Solves the adjusted system symmetric(matrixA) * solution = rhs using the MUMPS solver.
    mumps_solver_.solveInPlace(solution);
}

#endif
