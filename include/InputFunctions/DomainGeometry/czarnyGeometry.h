#pragma once

#include <cmath>

#include "../domainGeometry.h"

/* Triangular Shape and Ellipticity */

class CzarnyGeometry
{
public:
    explicit CzarnyGeometry();
    explicit CzarnyGeometry(double Rmax, double inverse_aspect_ratio_epsilon, double ellipticity_e);

    double Fx(double r, double theta) const;
    double Fy(double r, double theta) const;
    double dFx_dr(double r, double theta) const;
    double dFy_dr(double r, double theta) const;
    double dFx_dt(double r, double theta) const;
    double dFy_dt(double r, double theta) const;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};

#include "czarnyGeometry.inl"
