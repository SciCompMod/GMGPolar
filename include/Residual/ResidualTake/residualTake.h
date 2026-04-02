#pragma once

#include "../residual.h"

template <concepts::DomainGeometry DomainGeometry>
class ResidualTake : public Residual<DomainGeometry>
{
public:
    explicit ResidualTake(const PolarGrid& grid, const LevelCache<DomainGeometry>& level_cache,
                          const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior,
                          const int num_omp_threads);
    ~ResidualTake() override = default;

    void computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const override;

    void applySystemOperator(Vector<double> result, ConstVector<double> x) const override;

private:
    void applyCircleSection(const int i_r, Vector<double> result, ConstVector<double> x) const;
    void applyRadialSection(const int i_theta, Vector<double> result, ConstVector<double> x) const;
};

#include "residualTake.inl"
#include "applyATake.inl"
