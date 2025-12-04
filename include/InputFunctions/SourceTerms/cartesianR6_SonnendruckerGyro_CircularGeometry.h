#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR6_SonnendruckerGyro_CircularGeometry : public SourceTerm
{
public:
    CartesianR6_SonnendruckerGyro_CircularGeometry() = default;
    explicit CartesianR6_SonnendruckerGyro_CircularGeometry(const double& Rmax);
    virtual ~CartesianR6_SonnendruckerGyro_CircularGeometry() = default;

    double rhs_f(const double& r, const double& theta) const override;

private:
    const double Rmax = 1.3;
};
