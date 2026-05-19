#pragma once

template <class LevelCacheType>
ResidualTake<LevelCacheType>::ResidualTake(const PolarGrid<DefaultMemorySpace>& grid, const LevelCacheType& level_cache,
                                           bool DirBC_Interior, int num_omp_threads)
    : Residual<LevelCacheType>(grid, level_cache, DirBC_Interior, num_omp_threads)
{
}

/* ------------------ */
/* result = rhs - A*x */
template <class LevelCacheType>
void ResidualTake<LevelCacheType>::computeResidual(HostVector<double> h_result, HostConstVector<double> h_rhs,
                                                   HostConstVector<double> h_x) const
{
    assert(h_result.size() == h_x.size());

    applySystemOperator(h_result, h_x);

	auto rhs = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_rhs);
	auto result = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_result);

    // Subtract A*x from rhs to get the residual.
    const int n = result.size();

    Kokkos::parallel_for(
        "Residual Take: Subtract A*x from rhs", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
        KOKKOS_LAMBDA(const int i) { result[i] = rhs[i] - result[i]; });

    Kokkos::fence();

	Kokkos::deep_copy(h_result, result);
}
