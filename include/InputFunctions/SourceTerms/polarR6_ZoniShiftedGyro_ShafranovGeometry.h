#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_ZoniShiftedGyro_ShafranovGeometry : public SourceTerm
{
public:

    explicit PolarR6_ZoniShiftedGyro_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);
    virtual ~PolarR6_ZoniShiftedGyro_ShafranovGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
