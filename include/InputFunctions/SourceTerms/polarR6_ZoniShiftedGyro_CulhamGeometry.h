#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_ZoniShiftedGyro_CulhamGeometry : public SourceTerm
{
public:
    PolarR6_ZoniShiftedGyro_CulhamGeometry() = default;
    explicit PolarR6_ZoniShiftedGyro_CulhamGeometry(double Rmax);
    virtual ~PolarR6_ZoniShiftedGyro_CulhamGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
