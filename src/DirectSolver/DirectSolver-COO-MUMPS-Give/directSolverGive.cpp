#include "../../../include/DirectSolver/DirectSolver-COO-MUMPS-Give/directSolverGive.h"

#ifdef GMGPOLAR_USE_MUMPS

DirectSolver_COO_MUMPS_Give::DirectSolver_COO_MUMPS_Give(const PolarGrid& grid, const LevelCache& level_cache,
                                                         const DomainGeometry& domain_geometry,
                                                         const DensityProfileCoefficients& density_profile_coefficients,
                                                         bool DirBC_Interior, int num_omp_threads)
    : DirectSolver(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads)
{
    SparseMatrixCOO<double> solver_matrix = buildSolverMatrix();
    mumps_solver_.emplace(std::move(solver_matrix));
}

void DirectSolver_COO_MUMPS_Give::solveInPlace(Vector<double> solution)
{
    // Adjusts the right-hand side vector to account for symmetry corrections.
    // This transforms the system matrixA * solution = rhs into the equivalent system:
    // symmetric_DBc(matrixA) * solution = rhs - applySymmetryShift(rhs).
    // The correction modifies the rhs to account for the influence of the Dirichlet boundary conditions,
    // ensuring that the solution at the boundary is correctly adjusted and maintains the required symmetry.
    applySymmetryShift(solution);
    // Solves the adjusted system symmetric(matrixA) * solution = rhs using the MUMPS solver.
    mumps_solver_->solve(solution);
}

#endif
