#pragma once

#include "../residual.h"

namespace gmgpolar
{

template <class LevelCacheType>
class ResidualTake : public Residual<LevelCacheType>
{
public:
    explicit ResidualTake(const PolarGrid& grid, const LevelCacheType& level_cache, const bool DirBC_Interior,
                          const int num_omp_threads);
    KOKKOS_DEFAULTED_FUNCTION ResidualTake(const ResidualTake&) = default;
    KOKKOS_DEFAULTED_FUNCTION ~ResidualTake() override = default;

    void applySystemOperator(Vector<double> result, ConstVector<double> x) const final;
    void computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const final;
};

#include "residualTake.inl"
#include "applyATake.inl"

} // namespace gmgpolar
