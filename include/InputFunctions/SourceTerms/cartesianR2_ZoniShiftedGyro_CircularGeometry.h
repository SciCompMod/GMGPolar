#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR2_ZoniShiftedGyro_CircularGeometry : public SourceTerm
{
public:

    explicit CartesianR2_ZoniShiftedGyro_CircularGeometry(double Rmax);
    virtual ~CartesianR2_ZoniShiftedGyro_CircularGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
