#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR6_ZoniShifted_CircularGeometry : public SourceTerm { 
public:
    CartesianR6_ZoniShifted_CircularGeometry() = default;
    explicit CartesianR6_ZoniShifted_CircularGeometry(const double& Rmax);
    virtual ~CartesianR6_ZoniShifted_CircularGeometry() = default;

    double rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};
