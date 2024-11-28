#pragma once

#include "../extrapolatedSmoother.h"

class ExtrapolatedSmootherGive : public ExtrapolatedSmoother
{
public:
    explicit ExtrapolatedSmootherGive(const PolarGrid& grid,
                                      const LevelCache& level_cache,
                                      const DomainGeometry& domain_geometry,
                                      const DensityProfileCoefficients& density_profile_coefficients,
                                      bool DirBC_Interior,
                                      int num_omp_threads);

    ~ExtrapolatedSmootherGive() override;

    void extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) override;

private:
    void extrapolatedSmoothingInPlaceSequential(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void extrapolatedSmoothingInPlaceForLoop(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp); /* This is the fastest option */
    void extrapolatedSmoothingInPlaceTaskLoop(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void extrapolatedSmoothingInPlaceTaskDependencies(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);

    SparseMatrix<double> inner_boundary_circle_matrix_;
    DMUMPS_STRUC_C inner_boundary_mumps_solver_;
    std::vector<DiagonalSolver<double>> circle_diagonal_solver_;
    std::vector<DiagonalSolver<double>> radial_diagonal_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> circle_tridiagonal_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> radial_tridiagonal_solver_;

    // clang-format off
    Stencil stencil_center_ = 
        {-1, -1, -1,
        -1,  0, -1,
        -1, -1, -1};

    Stencil stencil_center_left_ = 
        {-1, -1, -1,
        1,  0, -1,
        -1, -1, -1};
    // clang-format on

    const Stencil& getStencil(int i_r, int i_theta) const;
    int getNonZeroCountCircleAsc(const int i_r) const;
    int getNonZeroCountRadialAsc(const int i_theta) const;

    int getCircleAscIndex(const int i_r, const int i_theta) const;
    int getRadialAscIndex(const int i_r, const int i_theta) const;

    void buildAscMatrices();
    void buildAscCircleSection(const int i_r);
    void buildAscRadialSection(const int i_theta);

    void applyAscOrthoCircleSection(const int i_r, const SmootherColor smoother_color, const Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color, const Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);

    void solveCircleSection(const int i_r, Vector<double>& x, Vector<double>& temp, Vector<double>& solver_storage_1, Vector<double>& solver_storage_2);
    void solveRadialSection(const int i_theta, Vector<double>& x, Vector<double>& temp, Vector<double>& solver_storage);

    void initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, const SparseMatrix<double>& solver_matrix);
    void finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver);
};
