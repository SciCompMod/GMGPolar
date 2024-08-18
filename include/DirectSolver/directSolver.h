#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <vector>
#include <iostream>

#include "mpi.h" 
#include "dmumps_c.h"   

#include "../PolarGrid/polargrid.h"
#include "../InputFunctions/domainGeometry.h"
#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/vector_operations.h"
#include "../common/constants.h"
#include "../Level/level.h"
#include "../Stencil/stencil.h"

class DirectSolver {
public:
    explicit DirectSolver(const PolarGrid& grid, const LevelCache& level_cache, 
                          const DomainGeometry& domain_geometry,
                          bool DirBC_Interior, int num_omp_threads);
    ~DirectSolver();

    void solveInPlace(Vector<double>& solution);
private:
    // Constructor parameters
    const PolarGrid& grid_;
    const std::vector<double>& sin_theta_cache_;
    const std::vector<double>& cos_theta_cache_;
    const std::vector<double>& coeff_alpha_cache_;
    const std::vector<double>& coeff_beta_cache_;
    const DomainGeometry& domain_geometry_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;

    // Solver matrix and MUMPS solver structure
    SparseMatrix<double> solver_matrix_;
    DMUMPS_STRUC_C mumps_solver_;

    // Constructs a symmetric solver matrix.
    SparseMatrix<double> buildSolverMatrix();
    void buildSolverMatrixCircleSection(const int i_r, SparseMatrix<double>& solver_matrix);
    void buildSolverMatrixRadialSection(const int i_theta, SparseMatrix<double>& solver_matrix);
    
    // Initializes the MUMPS solver with the specified matrix.
    void initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, const SparseMatrix<double>& solver_matrix);

    // Adjusts the right-hand side vector to account for symmetry corrections.
    // This transforms the system matrixA * solution = rhs into the equivalent system:
    // symmetric(A) * solution = rhs - applySymmetryShift(rhs).
    void applySymmetryShift(Vector<double>& rhs) const;
    void applySymmetryShiftInnerBoundary(Vector<double>& x) const;
    void applySymmetryShiftOuterBoundary(Vector<double>& x) const;
    
    // Solves the adjusted system symmetric(matrixA) * solution = rhs using the MUMPS solver.
    void solveWithMumps(Vector<double>& solution);

    // Finalizes the MUMPS solver, releasing any allocated resources.
    void finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver);

    // Returns the total number of non-zero elements in the solver matrix.
    int getNonZeroCountSolverMatrix() const;

    // Returns the index of the first non-zero element in the solver matrix for the given position.
    int getSolverMatrixIndex(int i_r, int i_theta) const;

    // Retrieves the stencil for the solver matrix at the given radial index.
    const Stencil& getStencil(int i_r) const;

};