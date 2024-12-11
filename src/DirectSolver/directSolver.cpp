#include "../../include/DirectSolver/directSolver.h"

DirectSolver::DirectSolver(const Level& level, const DomainGeometry& domain_geometry, 
    const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior)
    : grid_(level.grid())
    , level_cache_(level.levelCache())
    , domain_geometry_(domain_geometry)
    , density_profile_coefficients_(density_profile_coefficients)
    , DirBC_Interior_(DirBC_Interior)
{
    solver_matrix_ = buildSolverMatrix();
    initializeMumpsSolver(mumps_solver_, solver_matrix_);
}

void DirectSolver::solveInPlace(Vector<double>& solution)
{
    // Adjusts the right-hand side vector to account for symmetry corrections.
    // This transforms the system matrixA * solution = rhs into the equivalent system:
    // symmetric(matrixA) * solution = rhs - applySymmetryShift(rhs).
    applySymmetryShift(solution);
    // Solves the adjusted system symmetric(matrixA) * solution = rhs using the MUMPS solver.
    solveWithMumps(solution);
}

DirectSolver::~DirectSolver()
{
    finalizeMumpsSolver(mumps_solver_);
}