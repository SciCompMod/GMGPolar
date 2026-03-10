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

/* ------------ */
/* result = A*x */
template <concepts::DomainGeometry DomainGeometry>
void ResidualTake<DomainGeometry>::applySystemOperator(Vector<double> result, ConstVector<double> x) const
{
    assert(result.size() == x.size());

    const int num_omp_threads = Residual<DomainGeometry>::num_omp_threads_;

    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const PolarGrid& grid = Residual<DomainGeometry>::grid_;

#pragma omp parallel num_threads(num_omp_threads)
    {
/* Circle Section */
#pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            applyCircleSection(i_r, result, x);
        }

/* Radial Section */
#pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            applyRadialSection(i_theta, result, x);
        }
    }
}

/* ------------------ */
/* result = rhs - A*x */
template <concepts::DomainGeometry DomainGeometry>
void ResidualTake<DomainGeometry>::computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const
{
    assert(result.size() == x.size());

    const int num_omp_threads = Residual<DomainGeometry>::num_omp_threads_;

    applySystemOperator(result, x);

    // Subtract A*x from rhs to get the residual.
    const int n = result.size();
#pragma omp parallel for num_threads(num_omp_threads)
    for (int i = 0; i < n; i++) {
        result[i] = rhs[i] - result[i];
    }
}
