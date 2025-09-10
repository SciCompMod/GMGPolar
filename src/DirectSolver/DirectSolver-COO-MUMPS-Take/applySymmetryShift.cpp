#include "../../../include/DirectSolver/DirectSolver-COO-MUMPS-Take/directSolverTake.h"

#ifdef GMGPOLAR_USE_MUMPS

/* ----------------------- */
/* Boundary Symmetry Shift */
/* ----------------------- */

void DirectSolverTake::applySymmetryShiftInnerBoundary(
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x) const
{
    assert(DirBC_Interior_);

    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr = level_cache_.arr();
    const auto& att = level_cache_.att();
    const auto& art = level_cache_.art();

    const int i_r = 1;

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        double h1 = grid_.radialSpacing(i_r - 1);
        double k1 = grid_.angularSpacing(i_theta - 1);
        double k2 = grid_.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;

        const int i_theta_M1 = grid_.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid_.wrapThetaIndex(i_theta + 1);

        const int bottom_left = grid_.index(i_r - 1, i_theta_M1);
        const int left        = grid_.index(i_r - 1, i_theta);
        const int top_left    = grid_.index(i_r - 1, i_theta_P1);
        const int bottom      = grid_.index(i_r, i_theta_M1);
        const int center      = grid_.index(i_r, i_theta);
        const int top         = grid_.index(i_r, i_theta_P1);

        x[center] -= (-coeff1 * (arr[center] + arr[left]) * x[left] /* Left */
                      - 0.25 * (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */
                      + 0.25 * (art[left] + art[top]) * x[top_left] /* Top Left */
        );
    }
}

void DirectSolverTake::applySymmetryShiftOuterBoundary(
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x) const
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr = level_cache_.arr();
    const auto& att = level_cache_.att();
    const auto& art = level_cache_.art();

    const int i_r = grid_.nr() - 2;

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        double h2 = grid_.radialSpacing(i_r);
        double k1 = grid_.angularSpacing(i_theta - 1);
        double k2 = grid_.angularSpacing(i_theta);

        double coeff2 = 0.5 * (k1 + k2) / h2;

        const int i_theta_M1 = grid_.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid_.wrapThetaIndex(i_theta + 1);

        const int bottom       = grid_.index(i_r, i_theta_M1);
        const int center       = grid_.index(i_r, i_theta);
        const int top          = grid_.index(i_r, i_theta_P1);
        const int bottom_right = grid_.index(i_r + 1, i_theta_M1);
        const int right        = grid_.index(i_r + 1, i_theta);
        const int top_right    = grid_.index(i_r + 1, i_theta_P1);

        x[center] -= (-coeff2 * (arr[center] + arr[right]) * x[right] /* Right */
                      + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */
                      - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */
        );
    }
}

void DirectSolverTake::applySymmetryShift(Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x) const
{
    assert(x.size() == grid_.numberOfNodes());
    assert(grid_.nr() >= 4);

    if (num_omp_threads_ == 1) {
        /* Single-threaded execution */
        if (DirBC_Interior_) {
            applySymmetryShiftInnerBoundary(x);
        }
        applySymmetryShiftOuterBoundary(x);
    }
    else {
    #pragma omp parallel sections num_threads(num_omp_threads_)
        {
    #pragma omp section
            {
                if (DirBC_Interior_) {
                    applySymmetryShiftInnerBoundary(x);
                }
            }

    #pragma omp section
            {
                applySymmetryShiftOuterBoundary(x);
            }
        }
    }
}

#endif
