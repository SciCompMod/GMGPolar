#pragma once

#include <cmath>

#include "../domainGeometry.h"

class CircularGeometry
{
public:
    CircularGeometry() = default;
    explicit CircularGeometry(double Rmax);


    double Fx(double r, double theta) const;
    double Fy(double r, double theta) const;
    double dFx_dr(double r, double theta) const;
    double dFy_dr(double r, double theta) const;
    double dFx_dt(double r, double theta) const;
    double dFy_dt(double r, double theta) const;

private:
    const double Rmax = 1.3;
};

#include "circularGeometry.inl"