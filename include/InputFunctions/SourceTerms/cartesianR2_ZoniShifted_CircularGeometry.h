#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR2_ZoniShifted_CircularGeometry : public SourceTerm
{
public:

    explicit CartesianR2_ZoniShifted_CircularGeometry(double Rmax);
    virtual ~CartesianR2_ZoniShifted_CircularGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
