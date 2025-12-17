#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR6_ZoniShiftedGyro_CircularGeometry : public SourceTerm
{
public:

    explicit CartesianR6_ZoniShiftedGyro_CircularGeometry(double Rmax);
    virtual ~CartesianR6_ZoniShiftedGyro_CircularGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
