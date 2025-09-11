#pragma once

#include "../directSolver.h"
#include <Kokkos_Core.hpp>

class DirectSolverGiveCustomLU : public DirectSolver
{
public:
    explicit DirectSolverGiveCustomLU(const PolarGrid& grid, const LevelCache& level_cache,
                                      const DomainGeometry& domain_geometry,
                                      const DensityProfileCoefficients& density_profile_coefficients,
                                      bool DirBC_Interior, int num_omp_threads);

    ~DirectSolverGiveCustomLU() override;
    // Note: The rhs (right-hand side) vector gets overwritten with the solution.
    void solveInPlace(Vector<double> solution) override;

private:
    // Solver matrix and solver structure
    SparseMatrixCSR<double> solver_matrix_;
    SparseLUSolver<double> lu_solver_;

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
        7, 4, 8,
        1, 0, 2,
        5, 3, 6
    };
    const Stencil stencil_next_outer_DB_ = {
        7, 4, 8,
        1, 0, 2,
        5, 3, 6
    };
    // clang-format on

    SparseMatrixCSR<double> buildSolverMatrix();
    void buildSolverMatrixCircleSection(const int i_r, SparseMatrixCSR<double>& solver_matrix);
    void buildSolverMatrixRadialSection(const int i_theta, SparseMatrixCSR<double>& solver_matrix);

    // Returns the total number of non-zero elements in the solver matrix.
    int getNonZeroCountSolverMatrix() const;
    // Retrieves the stencil for the solver matrix at the given radial index.
    const Stencil& getStencil(int i_r) const;

    int getStencilSize(int global_index) const;
};
