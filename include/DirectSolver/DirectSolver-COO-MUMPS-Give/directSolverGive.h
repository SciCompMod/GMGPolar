#pragma once

#include "../directSolver.h"

#ifdef GMGPOLAR_USE_MUMPS

class DirectSolver_COO_MUMPS_Give : public DirectSolver
{
public:
    explicit DirectSolver_COO_MUMPS_Give(const PolarGrid& grid, const LevelCache& level_cache,
                                         const DomainGeometry& domain_geometry,
                                         const DensityProfileCoefficients& density_profile_coefficients,
                                         bool DirBC_Interior, int num_omp_threads);

    // Note: The rhs (right-hand side) vector gets overwritten during the solution process.
    void solveInPlace(Vector<double> solution) override;

private:
    // The stencil definitions must be defined before the declaration of the mumps_solver,
    // since the mumps solver will be build in the member initializer of the DirectSolver class.

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

    // MUMPS solver structure with the solver matrix initialized in the constructor.
    CooMumpsSolver mumps_solver_;

    // Constructs a symmetric solver matrix.
    SparseMatrixCOO<double> buildSolverMatrix();
    void buildSolverMatrixCircleSection(const int i_r, SparseMatrixCOO<double>& solver_matrix);
    void buildSolverMatrixRadialSection(const int i_theta, SparseMatrixCOO<double>& solver_matrix);

    // Adjusts the right-hand side vector for symmetry corrections.
    // This modifies the system from
    //    A * solution = rhs
    // to the equivalent system
    //    symmetric_DBc(A) * solution = rhs - applySymmetryShift(rhs).
    // The correction modifies the rhs to account for the influence of the Dirichlet boundary conditions,
    // ensuring that the solution at the boundary is correctly adjusted and maintains the required symmetry.
    void applySymmetryShift(Vector<double> rhs) const;
    void applySymmetryShiftInnerBoundary(Vector<double> x) const;
    void applySymmetryShiftOuterBoundary(Vector<double> x) const;

    // Returns the total number of non-zero elements in the solver matrix.
    int getNonZeroCountSolverMatrix() const;

    // Returns the index of the first non-zero element in the solver matrix for the given position.
    int getSolverMatrixIndex(int i_r, int i_theta) const;

    // Retrieves the stencil for the solver matrix at the given radial index.
    const Stencil& getStencil(int i_r) const;

    void nodeBuildSolverMatrixGive(int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                   SparseMatrixCOO<double>& solver_matrix, double arr, double att, double art,
                                   double detDF, double coeff_beta);
};

#endif
