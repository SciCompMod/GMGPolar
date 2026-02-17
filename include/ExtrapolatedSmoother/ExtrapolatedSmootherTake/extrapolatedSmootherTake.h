#pragma once

#include "../extrapolatedSmoother.h"

class ExtrapolatedSmootherTake : public ExtrapolatedSmoother
{
public:
    explicit ExtrapolatedSmootherTake(const PolarGrid& grid, const LevelCache& level_cache,
                                      const DomainGeometry& domain_geometry,
                                      const DensityProfileCoefficients& density_profile_coefficients,
                                      bool DirBC_Interior, int num_omp_threads);

    ~ExtrapolatedSmootherTake() override;

    void extrapolatedSmoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp) override;

private:
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
    using MatrixType = SparseMatrixCOO<double>;
    DMUMPS_STRUC_C inner_boundary_mumps_solver_;
#else
    using MatrixType = SparseMatrixCSR<double>;
    SparseLUSolver<double> inner_boundary_lu_solver_;
#endif
    MatrixType inner_boundary_circle_matrix_;

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
                                    ConstVector<double> rhs, Vector<double> temp);
    void applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color, ConstVector<double> x,
                                    ConstVector<double> rhs, Vector<double> temp);

    void solveCircleSection(const int i_r, Vector<double> x, Vector<double> temp, Vector<double> solver_storage_1,
                            Vector<double> solver_storage_2);
    void solveRadialSection(const int i_theta, Vector<double> x, Vector<double> temp, Vector<double> solver_storage);

#ifdef GMGPOLAR_USE_MUMPS
    void initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, SparseMatrixCOO<double>& solver_matrix);
    void finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver);
#endif

    void nodeBuildSmootherTake(int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                               MatrixType& inner_boundary_circle_matrix,
                               std::vector<DiagonalSolver<double>>& circle_diagonal_solver,
                               std::vector<DiagonalSolver<double>>& radial_diagonal_solver,
                               std::vector<SymmetricTridiagonalSolver<double>>& circle_tridiagonal_solver,
                               std::vector<SymmetricTridiagonalSolver<double>>& radial_tridiagonal_solver,
                               ConstVector<double>& arr, ConstVector<double>& att, ConstVector<double>& art,
                               ConstVector<double>& detDF, ConstVector<double>& coeff_beta);
};
