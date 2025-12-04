#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_ZoniShifted_CzarnyGeometry : public SourceTerm
{
public:
    PolarR6_ZoniShifted_CzarnyGeometry() = default;
    explicit PolarR6_ZoniShifted_CzarnyGeometry(double Rmax, double inverse_aspect_ratio_epsilon, double ellipticity_e);
    virtual ~PolarR6_ZoniShifted_CzarnyGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
