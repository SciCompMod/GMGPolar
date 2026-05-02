#pragma once

template <class LevelCacheType>
ResidualGive<LevelCacheType>::ResidualGive(const PolarGrid& grid, const LevelCacheType& level_cache,
                                           bool DirBC_Interior, int num_omp_threads)
    : Residual<LevelCacheType>(grid, level_cache, DirBC_Interior, num_omp_threads)
{
}

/* ------------------ */
/* result = rhs - A*x */
template <class LevelCacheType>
void ResidualGive<LevelCacheType>::computeResidual(Vector<double> result, ConstVector<double> rhs,
                                                   ConstVector<double> x) const
{
    assert(result.size() == x.size());

    applySystemOperator(result, x);

    // Subtract A*x from rhs to get the residual.
    const int n               = result.size();
    const int num_omp_threads = Residual<LevelCacheType>::num_omp_threads_;
#pragma omp parallel for num_threads(num_omp_threads)
    for (int i = 0; i < n; i++) {
        result[i] = rhs[i] - result[i];
    }
}
