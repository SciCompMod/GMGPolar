#pragma once

#include "../smoother.h"

#ifdef GMGPOLAR_USE_MUMPS
    #include "dmumps_c.h"
    #include "mpi.h"
#endif

class SmootherTake : public Smoother
{
public:
    explicit SmootherTake(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                          const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior,
                          int num_omp_threads);
    ~SmootherTake() override;

    void smoothing(Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x,
                   const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs,
                   Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> temp) override;

private:
    // The A_sc matrix on i_r = 0 is defined through the COO/CSR matrix
    // 'inner_boundary_circle_matrix_' due to the across-origin treatment.
    // It isn't tridiagonal and thus it requires a more advanced solver.
    // Note that circle_tridiagonal_solver_[0] is thus unused!
    // Additionally 'circle_tridiagonal_solver_[index]' will refer to the circular line i_r = index and
    // 'radial_tridiagonal_solver_[index] will refer to the radial line i_theta = index.
#ifdef GMGPOLAR_USE_MUMPS
    SparseMatrixCOO<double> inner_boundary_circle_matrix_;
    DMUMPS_STRUC_C inner_boundary_mumps_solver_;
#else
    SparseMatrixCSR<double> inner_boundary_circle_matrix_;
    SparseLUSolver<double> inner_boundary_lu_solver_;
#endif
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

    void applyAscOrthoCircleSection(const int i_r, const SmootherColor smoother_color,
                                    const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x,
                                    const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs,
                                    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> temp);
    void applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color,
                                    const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x,
                                    const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs,
                                    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> temp);

    void solveCircleSection(const int i_r, Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x,
                            Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> temp,
                            Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> solver_storage_1,
                            Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> solver_storage_2);
    void solveRadialSection(const int i_theta, Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x,
                            Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> temp,
                            Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> solver_storage);

#ifdef GMGPOLAR_USE_MUMPS
    void initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, SparseMatrixCOO<double>& solver_matrix);
    void finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver);
#endif
};
