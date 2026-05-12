#pragma once

#include <cmath>

#include "../domainGeometry.h"

/* Strechted Ellipse with a Shafranov Shift */

namespace gmgpolar
{

class ShafranovGeometry
{
public:
    ShafranovGeometry() = default;
    explicit ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);

    double Fx(double r, double theta) const;
    double Fy(double r, double theta) const;
    double dFx_dr(double r, double theta) const;
    double dFy_dr(double r, double theta) const;
    double dFx_dt(double r, double theta) const;
    double dFy_dt(double r, double theta) const;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};

#include "shafranovGeometry.inl"
} // namespace gmgpolar
