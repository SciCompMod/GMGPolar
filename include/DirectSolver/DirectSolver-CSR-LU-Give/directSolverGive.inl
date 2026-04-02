#pragma once

template <concepts::DomainGeometry DomainGeometry>
DirectSolver_CSR_LU_Give<DomainGeometry>::DirectSolver_CSR_LU_Give(const PolarGrid& grid,
                                                                   const LevelCache<DomainGeometry>& level_cache,
                                                                   bool DirBC_Interior, int num_omp_threads)
    : DirectSolver<DomainGeometry>(grid, level_cache, DirBC_Interior, num_omp_threads)
{
    solver_matrix_ = buildSolverMatrix();
    lu_solver_     = SparseLUSolver<double>(solver_matrix_);
}

template <concepts::DomainGeometry DomainGeometry>
void DirectSolver_CSR_LU_Give<DomainGeometry>::solveInPlace(Vector<double> solution)
{
    lu_solver_.solveInPlace(solution);
}
