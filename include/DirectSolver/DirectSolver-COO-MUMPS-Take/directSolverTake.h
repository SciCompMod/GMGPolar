#pragma once

#include "../directSolver.h"

#ifdef GMGPOLAR_USE_MUMPS

    #include "dmumps_c.h"
    #include "mpi.h"

class DirectSolverTake : public DirectSolver
{
public:
    explicit DirectSolverTake(const PolarGrid& grid, const LevelCache& level_cache,
                              const DomainGeometry& domain_geometry,
                              const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior,
                              int num_omp_threads);

    ~DirectSolverTake() override;
    // Note: The rhs (right-hand side) vector gets overwritten during the solution process.
    void solveInPlace(Vector<double> solution) override;

private:
    // Solver matrix and MUMPS solver structure
    SparseMatrixCOO<double> solver_matrix_;
    DMUMPS_STRUC_C mumps_solver_;

    // clang-format off
    const Stencil stencil_interior_      = {
        7, 4, 8,
        1, 0, 2,
        5, 3, 6
    };
    const Stencil stencil_across_origin_ = {
        -1, 4, 6,
        1, 0, 2,
        -1, 3, 5
    };
    const Stencil stencil_DB_            = {
        -1, -1, -1,
        -1,  0, -1,
        -1, -1, -1
    };
    const Stencil stencil_next_inner_DB_ = {
        -1,  3,  5,
        -1,  0,  1,
        -1,  2,  4
    };
    const Stencil stencil_next_outer_DB_ = {
        5,  3, -1,
        1,  0, -1,
        4,  2, -1
    };
    // clang-format on

    // Constructs a symmetric solver matrix.
    SparseMatrixCOO<double> buildSolverMatrix();
    void buildSolverMatrixCircleSection(const int i_r, SparseMatrixCOO<double>& solver_matrix);
    void buildSolverMatrixRadialSection(const int i_theta, SparseMatrixCOO<double>& solver_matrix);

    // Initializes the MUMPS solver with the specified matrix.
    // Converts to 1-based indexing.
    void initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, SparseMatrixCOO<double>& solver_matrix);

    // Adjusts the right-hand side vector for symmetry corrections.
    // This modifies the system from
    //    A * solution = rhs
    // to the equivalent system
    //    symmetric_DBc(A) * solution = rhs - applySymmetryShift(rhs).
    // The correction modifies the rhs to account for the influence of the Dirichlet boundary conditions,
    // ensuring that the solution at the boundary is correctly adjusted and maintains the required symmetry.
    void applySymmetryShift(Vector<double> const rhs) const;
    void applySymmetryShiftInnerBoundary(Vector<double> const x) const;
    void applySymmetryShiftOuterBoundary(Vector<double> const x) const;

    // Solves the adjusted system symmetric(matrixA) * solution = rhs using the MUMPS solver.
    void solveWithMumps(Vector<double> solution);

    // Finalizes the MUMPS solver, releasing any allocated resources.
    void finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver);

    // Returns the total number of non-zero elements in the solver matrix.
    int getNonZeroCountSolverMatrix() const;

    // Returns the index of the first non-zero element in the solver matrix for the given position.
    int getSolverMatrixIndex(int i_r, int i_theta) const;

    // Retrieves the stencil for the solver matrix at the given radial index.
    const Stencil& getStencil(int i_r) const;
};

#endif