#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR6_ZoniShifted_CircularGeometry : public SourceTerm
{
public:
    CartesianR6_ZoniShifted_CircularGeometry() = default;
    explicit CartesianR6_ZoniShifted_CircularGeometry(double Rmax);
    virtual ~CartesianR6_ZoniShifted_CircularGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
