#pragma once

#include "../extrapolatedSmoother.h"

#ifdef GMGPOLAR_USE_MUMPS
    #include "dmumps_c.h"
    #include "mpi.h"
#endif

class ExtrapolatedSmootherGive : public ExtrapolatedSmoother
{
public:
    explicit ExtrapolatedSmootherGive(const PolarGrid& grid, const LevelCache& level_cache,
                                      const DomainGeometry& domain_geometry,
                                      const DensityProfileCoefficients& density_profile_coefficients,
                                      bool DirBC_Interior, int num_omp_threads);

    ~ExtrapolatedSmootherGive() override;

    void extrapolatedSmoothing(Vector<double> const x, ConstVector<double> rhs, Vector<double> const temp) override;

private:
    void extrapolatedSmoothingSequential(Vector<double> const x, ConstVector<double> rhs, Vector<double> const temp);
    void extrapolatedSmoothingForLoop(Vector<double> const x, ConstVector<double> rhs, Vector<double> const temp);

    // The A_sc matrix on i_r = 0 is defined through the COO/CSR matrix
    // 'inner_boundary_circle_matrix_' due to the across-origin treatment.
    // It isn't tridiagonal and thus it requires a more advanced solver.
    // Note that circle_tridiagonal_solver_[0] is thus unused!

    // Lines containing coarse nodes are purely diagonal and thus are not stored in tridiagonal format.
    // - 'circle_tridiagonal_solver_[index] refers to the circular line i_r = 2*index,
    // - 'circle_diagonal_solver_[index] refers to the circular line i_r = 2*index + 1,
    // - 'radial_tridiagonal_solver_[index] refers to the radial line i_theta = 2*index,
    // - 'radial_diagonal_solver_[index] refers to the radial line i_theta = 2*index + 1.
#ifdef GMGPOLAR_USE_MUMPS
    SparseMatrixCOO<double> inner_boundary_circle_matrix_;
    DMUMPS_STRUC_C inner_boundary_mumps_solver_;
#else
    SparseMatrixCSR<double> inner_boundary_circle_matrix_;
    SparseLUSolver<double> inner_boundary_lu_solver_;
#endif
    std::vector<DiagonalSolver<double>> circle_diagonal_solver_;
    std::vector<DiagonalSolver<double>> radial_diagonal_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> circle_tridiagonal_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> radial_tridiagonal_solver_;

    // clang-format off
        Stencil stencil_center_ = {
        -1, -1, -1,
        -1,  0, -1,
        -1, -1, -1
    };
    Stencil stencil_center_left_ = {
        -1, -1, -1,
        1,  0, -1,
        -1, -1, -1
    };
    // clang-format on

    const Stencil& getStencil(int i_r, int i_theta) const;
    int getNonZeroCountCircleAsc(const int i_r) const;
    int getNonZeroCountRadialAsc(const int i_theta) const;

    int getCircleAscIndex(const int i_r, const int i_theta) const;
    int getRadialAscIndex(const int i_r, const int i_theta) const;

    void buildAscMatrices();
    void buildAscCircleSection(const int i_r);
    void buildAscRadialSection(const int i_theta);

    void applyAscOrthoCircleSection(const int i_r, const SmootherColor smoother_color, ConstVector<double> x,
                                    ConstVector<double> rhs, Vector<double> const temp);
    void applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color, ConstVector<double> x,
                                    ConstVector<double> rhs, Vector<double> const temp);

    void solveCircleSection(const int i_r, Vector<double> const x, Vector<double> const temp,
                            Vector<double> const solver_storage_1, Vector<double> const solver_storage_2);
    void solveRadialSection(const int i_theta, Vector<double> const x, Vector<double> const temp,
                            Vector<double> const solver_storage);

#ifdef GMGPOLAR_USE_MUMPS
    void initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, SparseMatrixCOO<double>& solver_matrix);
    void finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver);
#endif
};
