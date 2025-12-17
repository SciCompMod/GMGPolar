#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_ZoniGyro_CircularGeometry : public SourceTerm
{
public:

    explicit PolarR6_ZoniGyro_CircularGeometry(double Rmax);
    virtual ~PolarR6_ZoniGyro_CircularGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
