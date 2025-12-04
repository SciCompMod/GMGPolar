#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR6_SonnendruckerGyro_CzarnyGeometry : public SourceTerm
{
public:
    CartesianR6_SonnendruckerGyro_CzarnyGeometry() = default;
    explicit CartesianR6_SonnendruckerGyro_CzarnyGeometry(double Rmax, double inverse_aspect_ratio_epsilon,
                                                          double ellipticity_e);
    virtual ~CartesianR6_SonnendruckerGyro_CzarnyGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
