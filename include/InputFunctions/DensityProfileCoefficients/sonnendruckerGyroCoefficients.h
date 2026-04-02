#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class SonnendruckerGyroCoefficients
{
public:
    SonnendruckerGyroCoefficients() = default;
    explicit SonnendruckerGyroCoefficients(double Rmax, double alpha);

    double alpha(double r, double theta) const;
    double beta(double r, double theta) const;

    double getAlphaJump() const;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.66 * 1.3;
};
