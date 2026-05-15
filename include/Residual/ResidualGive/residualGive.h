#pragma once

#include "../residual.h"

namespace gmgpolar
{

template <class LevelCacheType>
class ResidualGive : public Residual<LevelCacheType>
{
public:
    explicit ResidualGive(const PolarGrid& grid, const LevelCacheType& level_cache, const bool DirBC_Interior);
    ~ResidualGive() override = default;

    void applySystemOperator(Vector<double> result, ConstVector<double> x) const final;
    void computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const final;
};

#include "residualGive.inl"
#include "applyAGive.inl"

} // namespace gmgpolar
