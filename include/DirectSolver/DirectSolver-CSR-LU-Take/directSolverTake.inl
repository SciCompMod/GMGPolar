#pragma once

template <concepts::DomainGeometry DomainGeometry>
DirectSolver_CSR_LU_Take<DomainGeometry>::DirectSolver_CSR_LU_Take(
    const PolarGrid& grid, const LevelCache<DomainGeometry>& level_cache, const DomainGeometry& domain_geometry,
    const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior, int num_omp_threads)
    : DirectSolver<DomainGeometry>(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior,
                                   num_omp_threads)
{
    solver_matrix_ = buildSolverMatrix();
    lu_solver_     = SparseLUSolver<double>(solver_matrix_);
}

template <concepts::DomainGeometry DomainGeometry>
void DirectSolver_CSR_LU_Take<DomainGeometry>::solveInPlace(Vector<double> solution)
{
    // Adjusts the right-hand side vector to account for symmetry corrections.
    // This transforms the system matrixA * solution = rhs into the equivalent system:
    // symmetric_DBc(matrixA) * solution = rhs - applySymmetryShift(rhs).
    // The correction modifies the rhs to account for the influence of the Dirichlet boundary conditions,
    // ensuring that the solution at the boundary is correctly adjusted and maintains the required symmetry.
    applySymmetryShift(solution);
    // Solves the adjusted system symmetric(matrixA) * solution = rhs using the MUMPS solver.
    lu_solver_.solveInPlace(solution);
}
