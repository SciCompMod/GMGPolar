#pragma once

#include "../residual.h"

class ResidualGive : public Residual
{
public:
    explicit ResidualGive(const PolarGrid& grid,
                          const LevelCache& level_cache,
                          const DomainGeometry& domain_geometry,
                          const DensityProfileCoefficients& density_profile_coefficients,
                          const bool DirBC_Interior,
                          const int num_omp_threads);
    ~ResidualGive() override = default;

    void computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const override;

private:
    void applyCircleSection(const int i_r, Vector<double>& result, const Vector<double>& x, const double& factor) const;
    void applyRadialSection(const int i_theta, Vector<double>& result, const Vector<double>& x, const double& factor) const;
};