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
    ~ResidualTake() override = default;

    void computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const override;

    void applySystemOperator(Vector<double> result, ConstVector<double> x) const override;
};

#include "residualTake.inl"
#include "applyATake.inl"
} // namespace gmgpolar
