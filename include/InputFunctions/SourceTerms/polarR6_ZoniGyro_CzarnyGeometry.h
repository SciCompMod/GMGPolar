#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_ZoniGyro_CzarnyGeometry : public SourceTerm
{
public:
    PolarR6_ZoniGyro_CzarnyGeometry() = default;
    explicit PolarR6_ZoniGyro_CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon,
                                             const double& ellipticity_e);
    virtual ~PolarR6_ZoniGyro_CzarnyGeometry() = default;

    double rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
