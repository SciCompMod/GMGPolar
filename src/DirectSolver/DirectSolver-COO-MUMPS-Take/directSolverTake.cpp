#include "../../../include/DirectSolver/DirectSolver-COO-MUMPS-Take/directSolverTake.h"

#ifdef GMGPOLAR_USE_MUMPS

DirectSolverTake::DirectSolverTake(const PolarGrid& grid, const LevelCache& level_cache,
                                   const DomainGeometry& domain_geometry,
                                   const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior,
                                   int num_omp_threads)
    : DirectSolver(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads)
{
    solver_matrix_ = buildSolverMatrix();
    initializeMumpsSolver(mumps_solver_, solver_matrix_);
}

void DirectSolverTake::solveInPlace(Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> solution)
{
    // Adjusts the right-hand side vector to account for symmetry corrections.
    // This transforms the system matrixA * solution = rhs into the equivalent system:
    // symmetric_DBc(matrixA) * solution = rhs - applySymmetryShift(rhs).
    // The correction modifies the rhs to account for the influence of the Dirichlet boundary conditions,
    // ensuring that the solution at the boundary is correctly adjusted and maintains the required symmetry.
    applySymmetryShift(solution);
    // Solves the adjusted system symmetric(matrixA) * solution = rhs using the MUMPS solver.
    solveWithMumps(solution);
}

DirectSolverTake::~DirectSolverTake()
{
    finalizeMumpsSolver(mumps_solver_);
}

#endif
