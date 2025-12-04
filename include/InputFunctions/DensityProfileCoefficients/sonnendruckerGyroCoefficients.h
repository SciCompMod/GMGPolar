#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class SonnendruckerGyroCoefficients : public DensityProfileCoefficients
{
public:
    SonnendruckerGyroCoefficients() = default;
    explicit SonnendruckerGyroCoefficients(double Rmax, double alpha);
    virtual ~SonnendruckerGyroCoefficients() = default;

    double alpha(double r, double theta) const override;
    double beta(double r, double theta) const override;

    double getAlphaJump() const override;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.66 * 1.3;
};
