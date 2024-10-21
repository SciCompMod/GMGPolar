#pragma once

#include "../directSolver.h"

class DirectSolverGive : public DirectSolver {
public:
    explicit DirectSolverGive(
        const PolarGrid& grid, const LevelCache& level_cache, 
        const DomainGeometry& domain_geometry, const DensityProfileCoefficients& density_profile_coefficients,
        bool DirBC_Interior, int num_omp_threads
    );
    ~DirectSolverGive() override;

    void solveInPlace(Vector<double>& solution) override;

private:
    // Solver matrix and MUMPS solver structure
    SparseMatrix<double> solver_matrix_;
    DMUMPS_STRUC_C mumps_solver_;

    const Stencil stencil_interior_ = 
        {7, 4, 8,
        1, 0, 2,
        5, 3, 6};
    const Stencil stencil_across_origin_ = 
        {-1, 4, 6,
        1, 0, 2,
        -1, 3, 5};
    const Stencil stencil_DB_ = 
        {-1, -1, -1,
        -1,  0, -1,
        -1, -1, -1};
    const Stencil stencil_next_inner_DB_ = 
        {-1, 3, 5,
        -1, 0, 1,
        -1, 2, 4};
    const Stencil stencil_next_outer_DB_ = 
        {5, 3, -1,
        1, 0, -1,
        4, 2, -1};

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
