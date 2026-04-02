#pragma once

template <class LevelCacheType>
DirectSolver_CSR_LU_Take<LevelCacheType>::DirectSolver_CSR_LU_Take(const PolarGrid& grid,
                                                                   const LevelCacheType& level_cache,
                                                                   bool DirBC_Interior, int num_omp_threads)
    : DirectSolver<LevelCacheType>(grid, level_cache, DirBC_Interior, num_omp_threads)
{
    solver_matrix_ = buildSolverMatrix();
    lu_solver_     = SparseLUSolver<double>(solver_matrix_);
}

template <class LevelCacheType>
void DirectSolver_CSR_LU_Take<LevelCacheType>::solveInPlace(Vector<double> solution)
{
    lu_solver_.solveInPlace(solution);
}
