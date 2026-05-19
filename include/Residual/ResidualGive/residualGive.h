#pragma once

#include "../residual.h"

namespace gmgpolar
{

template <class LevelCacheType>
class ResidualGive : public Residual<LevelCacheType>
{
public:
    explicit ResidualGive(const PolarGrid<Kokkos::HostSpace>& grid, const LevelCacheType& level_cache, const bool DirBC_Interior,
                          const int num_omp_threads);
    ~ResidualGive() override = default;

    void applySystemOperator(HostVector<double> result, HostConstVector<double> x) const final;
    void computeResidual(HostVector<double> result, HostConstVector<double> rhs, HostConstVector<double> x) const final;
};

#include "residualGive.inl"
#include "applyAGive.inl"

} // namespace gmgpolar
