#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_ZoniShiftedGyro_CulhamGeometry : public SourceTerm
{
public:
    PolarR6_ZoniShiftedGyro_CulhamGeometry() = default;
    explicit PolarR6_ZoniShiftedGyro_CulhamGeometry(const double& Rmax);
    virtual ~PolarR6_ZoniShiftedGyro_CulhamGeometry() = default;

    double rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};
