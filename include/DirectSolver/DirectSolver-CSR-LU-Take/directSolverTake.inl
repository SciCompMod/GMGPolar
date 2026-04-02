#pragma once

template <concepts::DomainGeometry DomainGeometry>
DirectSolver_CSR_LU_Take<DomainGeometry>::DirectSolver_CSR_LU_Take(
    const PolarGrid& grid, const LevelCache<DomainGeometry>& level_cache,
    const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior, int num_omp_threads)
    : DirectSolver<DomainGeometry>(grid, level_cache, density_profile_coefficients, DirBC_Interior, num_omp_threads)
{
    solver_matrix_ = buildSolverMatrix();
    lu_solver_     = SparseLUSolver<double>(solver_matrix_);
}

template <concepts::DomainGeometry DomainGeometry>
void DirectSolver_CSR_LU_Take<DomainGeometry>::solveInPlace(Vector<double> solution)
{
    lu_solver_.solveInPlace(solution);
}
