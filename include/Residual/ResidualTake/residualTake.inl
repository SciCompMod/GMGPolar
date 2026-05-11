#pragma once

template <class LevelCacheType>
ResidualTake<LevelCacheType>::ResidualTake(const PolarGrid& grid, const LevelCacheType& level_cache,
                                           bool DirBC_Interior, int num_omp_threads)
    : Residual<LevelCacheType>(grid, level_cache, DirBC_Interior, num_omp_threads)
{
}

/* ------------------ */
/* result = rhs - A*x */
template <class LevelCacheType>
void ResidualTake<LevelCacheType>::computeResidual(Vector<double> result, ConstVector<double> rhs,
                                                   ConstVector<double> x) const
{
    assert(result.size() == x.size());

    applySystemOperator(result, x);

    // Subtract A*x from rhs to get the residual.
    const int n = result.size();

    Kokkos::parallel_for(
        "Residual Take: Subtract A*x from rhs", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(const int i) { result[i] = rhs[i] - result[i]; });
}
