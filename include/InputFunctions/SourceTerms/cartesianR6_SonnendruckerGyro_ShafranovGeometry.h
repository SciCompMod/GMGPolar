#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR6_SonnendruckerGyro_ShafranovGeometry : public SourceTerm
{
public:

    explicit CartesianR6_SonnendruckerGyro_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);
    virtual ~CartesianR6_SonnendruckerGyro_ShafranovGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
