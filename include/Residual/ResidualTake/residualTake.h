#pragma once

#include "../residual.h"

namespace gmgpolar
{

template <class LevelCacheType>
class ResidualTake : public Residual<LevelCacheType>
{
public:
    explicit ResidualTake(const PolarGrid<Kokkos::HostSpace>& grid, const LevelCacheType& level_cache, const bool DirBC_Interior,
                          const int num_omp_threads);
    KOKKOS_DEFAULTED_FUNCTION ResidualTake(const ResidualTake&) = default;
    KOKKOS_DEFAULTED_FUNCTION ~ResidualTake() override          = default;

    void applySystemOperator(HostVector<double> result, HostConstVector<double> x) const final;
    void computeResidual(HostVector<double> result, HostConstVector<double> rhs, HostConstVector<double> x) const final;
};

#include "residualTake.inl"
#include "applyATake.inl"

} // namespace gmgpolar
