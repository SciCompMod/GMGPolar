#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class SonnendruckerCoefficients
{
public:
    SonnendruckerCoefficients() = default;
    explicit SonnendruckerCoefficients(double Rmax, double alpha);

    double alpha(double r, double theta) const;
    double beta(double r, double theta) const;

    double getAlphaJump() const;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.66 * 1.3;
};
