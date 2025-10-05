#pragma once

#include "../system_operator.h"

class SystemOperatorTake : public SystemOperator
{
public:
    explicit SystemOperatorTake(const PolarGrid& grid, const LevelCache& level_cache,
                                const DomainGeometry& domain_geometry,
                                const DensityProfileCoefficients& density_profile_coefficients,
                                const bool DirBC_Interior, const int num_omp_threads);
    ~SystemOperatorTake() override = default;

    void apply(Vector<double>& result, const Vector<double>& x) const override;

private:
    void applyCircleSection(const int i_r, Vector<double>& result, const Vector<double>& x) const;
    void applyRadialSection(const int i_theta, Vector<double>& result, const Vector<double>& x) const;
};
