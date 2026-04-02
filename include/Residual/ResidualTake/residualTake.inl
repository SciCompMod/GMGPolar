#pragma once

template <class LevelCacheType>
ResidualTake<LevelCacheType>::ResidualTake(const PolarGrid& grid, const LevelCacheType& level_cache,
                                           bool DirBC_Interior, int num_omp_threads)
    : Residual<LevelCacheType>(grid, level_cache, DirBC_Interior, num_omp_threads)
{
}

/* ------------ */
/* result = A*x */
template <class LevelCacheType>
void ResidualTake<LevelCacheType>::applySystemOperator(Vector<double> result, ConstVector<double> x) const
{
    assert(result.size() == x.size());

    const int num_omp_threads = Residual<LevelCacheType>::num_omp_threads_;

    assert(Residual<LevelCacheType>::level_cache_.cacheDensityProfileCoefficients());
    assert(Residual<LevelCacheType>::level_cache_.cacheDomainGeometry());

    const PolarGrid& grid = Residual<LevelCacheType>::grid_;

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
template <class LevelCacheType>
void ResidualTake<LevelCacheType>::computeResidual(Vector<double> result, ConstVector<double> rhs,
                                                   ConstVector<double> x) const
{
    assert(result.size() == x.size());

    const int num_omp_threads = Residual<LevelCacheType>::num_omp_threads_;

    applySystemOperator(result, x);

    // Subtract A*x from rhs to get the residual.
    const int n = result.size();
#pragma omp parallel for num_threads(num_omp_threads)
    for (int i = 0; i < n; i++) {
        result[i] = rhs[i] - result[i];
    }
}
