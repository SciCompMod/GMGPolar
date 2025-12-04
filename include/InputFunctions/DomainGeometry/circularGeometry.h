#pragma once

#include <cmath>

#include "../domainGeometry.h"

class CircularGeometry : public DomainGeometry
{
public:
    CircularGeometry() = default;
    explicit CircularGeometry(double Rmax);

    virtual ~CircularGeometry() = default;

    double Fx(double r, double theta) const override;
    double Fy(double r, double theta) const override;
    double dFx_dr(double r, double theta) const override;
    double dFy_dr(double r, double theta) const override;
    double dFx_dt(double r, double theta) const override;
    double dFy_dt(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};

#include "circularGeometry.inl"