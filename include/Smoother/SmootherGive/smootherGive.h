#pragma once

#include "../smoother.h"

class SmootherGive : public Smoother
{
public:
    explicit SmootherGive(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                          const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior,
                          int num_omp_threads);
    ~SmootherGive() override;

    void smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) override;

private:
    void smoothingInPlaceSequential(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void smoothingInPlaceForLoop(Vector<double>& x, const Vector<double>& rhs,
                                 Vector<double>& temp); /* This is the fastest option */
    void smoothingInPlaceTaskLoop(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void smoothingInPlaceTaskDependencies(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);

    SparseMatrixCOO<double> inner_boundary_circle_matrix_;
    DMUMPS_STRUC_C inner_boundary_mumps_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> circle_tridiagonal_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> radial_tridiagonal_solver_;

    // clang-format off
    const Stencil stencil_DB_ = {
        -1, -1, -1,
        -1,  0, -1,
        -1, -1, -1
    };
    /* Circle Stencils */
    const Stencil circle_stencil_interior_ = {
        -1,  2, -1,
        -1,  0, -1,
        -1,  1, -1
    };
    const Stencil circle_stencil_across_origin_ = {
        -1,  3, -1,
        1,  0, -1,
        -1,  2, -1
    };
    /* Radial Stencils */
    const Stencil radial_stencil_interior_ = {
        -1, -1, -1,
        1,  0,  2,
        -1, -1, -1
    };
    const Stencil radial_stencil_next_outer_DB_ = {
        -1, -1, -1,
        1,  0, -1,
        -1, -1, -1
    };
    const Stencil radial_stencil_next_circular_smoothing_ = {
        -1, -1, -1,
        -1,  0,  1,
        -1, -1, -1
    };
    // clang-format on

    const Stencil& getStencil(int i_r) const;
    int getNonZeroCountCircleAsc(const int i_r) const;
    int getNonZeroCountRadialAsc(const int i_theta) const;

    int getCircleAscIndex(const int i_r, const int i_theta) const;
    int getRadialAscIndex(const int i_r, const int i_theta) const;

    void buildAscMatrices();
    void buildAscCircleSection(const int i_r);
    void buildAscRadialSection(const int i_theta);

    void applyAscOrthoCircleSection(const int i_r, const SmootherColor smoother_color, const Vector<double>& x,
                                    const Vector<double>& rhs, Vector<double>& temp);
    void applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color, const Vector<double>& x,
                                    const Vector<double>& rhs, Vector<double>& temp);

    void solveCircleSection(const int i_r, Vector<double>& x, Vector<double>& temp, Vector<double>& solver_storage_1,
                            Vector<double>& solver_storage_2);
    void solveRadialSection(const int i_theta, Vector<double>& x, Vector<double>& temp, Vector<double>& solver_storage);

    void initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, const SparseMatrixCOO<double>& solver_matrix);
    void finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver);
};
