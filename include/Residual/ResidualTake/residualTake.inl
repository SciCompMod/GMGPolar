#pragma once

template <concepts::DomainGeometry DomainGeometry>
ResidualTake<DomainGeometry>::ResidualTake(const PolarGrid& grid, const LevelCache<DomainGeometry>& level_cache,
                                           const DomainGeometry& domain_geometry,
                                           const DensityProfileCoefficients& density_profile_coefficients,
                                           bool DirBC_Interior, int num_omp_threads)
    : Residual<DomainGeometry>(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior,
                               num_omp_threads)
{
}

/* ------------------ */
/* result = rhs - A*x */

// clang-format off
template <concepts::DomainGeometry DomainGeometry>
void ResidualTake<DomainGeometry>::computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const
{
    assert(result.size() == x.size());

    assert(Residual<DomainGeometry>::level_cache_.cacheDensityProfileCoefficients());
    assert(Residual<DomainGeometry>::level_cache_.cacheDomainGeometry());

    const PolarGrid& grid            = Residual<DomainGeometry>::grid_;
    const int        num_omp_threads = Residual<DomainGeometry>::num_omp_threads_;

    #pragma omp parallel num_threads(num_omp_threads)
    {
        /* Circle Section */
        #pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            applyCircleSection(i_r, result, rhs, x);
        }

        /* Radial Section */
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            applyRadialSection(i_theta, result, rhs, x);
        }
    }
}
// clang-format on
