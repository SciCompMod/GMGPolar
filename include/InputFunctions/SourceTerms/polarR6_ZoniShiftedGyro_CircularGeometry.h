#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_ZoniShiftedGyro_CircularGeometry : public SourceTerm
{
public:
    PolarR6_ZoniShiftedGyro_CircularGeometry() = default;
    explicit PolarR6_ZoniShiftedGyro_CircularGeometry(const double& Rmax);
    virtual ~PolarR6_ZoniShiftedGyro_CircularGeometry() = default;

    double rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};
