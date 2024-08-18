#include "../../include/DirectSolver/directSolver.h"

DirectSolver::DirectSolver(const PolarGrid& grid, const LevelCache& level_cache, 
                           const DomainGeometry& domain_geometry,
                           bool DirBC_Interior, int num_omp_threads
) :
    grid_(grid),
    sin_theta_cache_(level_cache.sin_theta()),
    cos_theta_cache_(level_cache.cos_theta()),
    coeff_alpha_cache_(level_cache.coeff_alpha()),
    coeff_beta_cache_(level_cache.coeff_beta()),
    domain_geometry_(domain_geometry),
    DirBC_Interior_(DirBC_Interior),
    num_omp_threads_(num_omp_threads)
{
    solver_matrix_ = buildSolverMatrix();
    initializeMumpsSolver(mumps_solver_, solver_matrix_);
}

void DirectSolver::solveInPlace(Vector<double>& solution) {
    // Adjusts the right-hand side vector to account for symmetry corrections.
    // This transforms the system matrixA * solution = rhs into the equivalent system:
    // symmetric(matrixA) * solution = rhs - applySymmetryShift(rhs).
    applySymmetryShift(solution);
    // Solves the adjusted system symmetric(matrixA) * solution = rhs using the MUMPS solver.
    solveWithMumps(solution);
}


DirectSolver::~DirectSolver() {
    finalizeMumpsSolver(mumps_solver_);
}