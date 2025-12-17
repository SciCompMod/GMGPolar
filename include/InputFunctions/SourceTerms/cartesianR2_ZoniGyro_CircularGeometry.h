#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR2_ZoniGyro_CircularGeometry : public SourceTerm
{
public:

    explicit CartesianR2_ZoniGyro_CircularGeometry(PolarGrid const& grid, double Rmax);
    virtual ~CartesianR2_ZoniGyro_CircularGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
