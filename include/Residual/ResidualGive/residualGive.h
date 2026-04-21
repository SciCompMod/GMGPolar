#pragma once

#include "../residual.h"

namespace gmgpolar
{

template <class LevelCacheType>
class ResidualGive : public Residual<LevelCacheType>
{
public:
    explicit ResidualGive(const PolarGrid& grid, const LevelCacheType& level_cache, const bool DirBC_Interior,
                          const int num_omp_threads);
    ~ResidualGive() override = default;

    void computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const final;

    void applySystemOperator(Vector<double> result, ConstVector<double> x) const final;

private:
    void applyCircleSection(const int i_r, Vector<double> result, ConstVector<double> x) const;
    void applyRadialSection(const int i_theta, Vector<double> result, ConstVector<double> x) const;
};

#include "residualGive.inl"
#include "applyAGive.inl"
} // namespace gmgpolar
