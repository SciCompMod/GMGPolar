#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_SonnendruckerGyro_CircularGeometry : public SourceTerm
{
public:
    PolarR6_SonnendruckerGyro_CircularGeometry() = default;
    explicit PolarR6_SonnendruckerGyro_CircularGeometry(double Rmax);
    virtual ~PolarR6_SonnendruckerGyro_CircularGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
