#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR2_SonnendruckerGyro_ShafranovGeometry : public SourceTerm
{
public:
    CartesianR2_SonnendruckerGyro_ShafranovGeometry() = default;
    explicit CartesianR2_SonnendruckerGyro_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);
    virtual ~CartesianR2_SonnendruckerGyro_ShafranovGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
