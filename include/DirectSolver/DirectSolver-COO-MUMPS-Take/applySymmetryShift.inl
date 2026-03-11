#pragma once

#ifdef GMGPOLAR_USE_MUMPS

/* ----------------------- */
/* Boundary Symmetry Shift */
/* ----------------------- */

template <concepts::DomainGeometry DomainGeometry>
void DirectSolver_COO_MUMPS_Take<DomainGeometry>::applySymmetryShiftInnerBoundary(Vector<double> x) const
{
    const PolarGrid& grid                         = DirectSolver<DomainGeometry>::grid_;
    const LevelCache<DomainGeometry>& level_cache = DirectSolver<DomainGeometry>::level_cache_;

    assert(DirectSolver<DomainGeometry>::DirBC_Interior_);

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    ConstVector<double> arr = level_cache.arr();
    ConstVector<double> att = level_cache.att();
    ConstVector<double> art = level_cache.art();

    const int i_r = 1;

    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
        double h1 = grid.radialSpacing(i_r - 1);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int bottom_left = grid.index(i_r - 1, i_theta_M1);
        const int left        = grid.index(i_r - 1, i_theta);
        const int top_left    = grid.index(i_r - 1, i_theta_P1);
        const int bottom      = grid.index(i_r, i_theta_M1);
        const int center      = grid.index(i_r, i_theta);
        const int top         = grid.index(i_r, i_theta_P1);

        x[center] -= (-coeff1 * (arr[center] + arr[left]) * x[left] /* Left */
                      - 0.25 * (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */
                      + 0.25 * (art[left] + art[top]) * x[top_left] /* Top Left */
        );
    }
}

template <concepts::DomainGeometry DomainGeometry>
void DirectSolver_COO_MUMPS_Take<DomainGeometry>::applySymmetryShiftOuterBoundary(Vector<double> x) const
{
    const PolarGrid& grid                         = DirectSolver<DomainGeometry>::grid_;
    const LevelCache<DomainGeometry>& level_cache = DirectSolver<DomainGeometry>::level_cache_;

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    ConstVector<double> arr = level_cache.arr();
    ConstVector<double> att = level_cache.att();
    ConstVector<double> art = level_cache.art();

    const int i_r = grid.nr() - 2;

    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff2 = 0.5 * (k1 + k2) / h2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int bottom       = grid.index(i_r, i_theta_M1);
        const int center       = grid.index(i_r, i_theta);
        const int top          = grid.index(i_r, i_theta_P1);
        const int bottom_right = grid.index(i_r + 1, i_theta_M1);
        const int right        = grid.index(i_r + 1, i_theta);
        const int top_right    = grid.index(i_r + 1, i_theta_P1);

        x[center] -= (-coeff2 * (arr[center] + arr[right]) * x[right] /* Right */
                      + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */
                      - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */
        );
    }
}

template <concepts::DomainGeometry DomainGeometry>
void DirectSolver_COO_MUMPS_Take<DomainGeometry>::applySymmetryShift(Vector<double> x) const
{
    const PolarGrid& grid     = DirectSolver<DomainGeometry>::grid_;
    const bool DirBC_Interior = DirectSolver<DomainGeometry>::DirBC_Interior_;
    const int num_omp_threads = DirectSolver<DomainGeometry>::num_omp_threads_;

    assert(std::ssize(x) == grid.numberOfNodes());
    assert(grid.nr() >= 4);

    #pragma omp parallel sections num_threads(num_omp_threads)
    {
    #pragma omp section
        {
            if (DirBC_Interior) {
                applySymmetryShiftInnerBoundary(x);
            }
        }

    #pragma omp section
        {
            applySymmetryShiftOuterBoundary(x);
        }
    }
}

#endif
