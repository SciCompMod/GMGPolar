#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR2_ZoniGyro_CircularGeometry : public SourceTerm { 
public:
    CartesianR2_ZoniGyro_CircularGeometry() = default;
    explicit CartesianR2_ZoniGyro_CircularGeometry(const double& Rmax);
    virtual ~CartesianR2_ZoniGyro_CircularGeometry() = default;

    double rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};
