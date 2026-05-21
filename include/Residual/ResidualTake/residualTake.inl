#pragma once

template <class LevelCacheType>
ResidualTake<LevelCacheType>::ResidualTake(const PolarGrid& grid, const LevelCacheType& level_cache,
                                           bool DirBC_Interior)
    : Residual<LevelCacheType>(grid, level_cache, DirBC_Interior)
{
}

/* ------------------ */
/* result = rhs - A*x */
template <class LevelCacheType>
void ResidualTake<LevelCacheType>::computeResidual(HostVector<double> result, HostConstVector<double> rhs,
                                                   HostConstVector<double> x) const
{
    assert(result.size() == x.size());

    applySystemOperator(result, x);

    // Subtract A*x from rhs to get the residual.
    const int n = result.size();

    Kokkos::parallel_for(
        "Residual Take: Subtract A*x from rhs", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
        KOKKOS_LAMBDA(const int i) { result[i] = rhs[i] - result[i]; });

    Kokkos::fence();
}
