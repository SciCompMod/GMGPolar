#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class PoissonCoefficients
{
public:
    PoissonCoefficients() = default;
    explicit PoissonCoefficients(double Rmax, double alpha);

    double alpha(double r, double theta) const;
    double beta(double r, double theta) const;

    double getAlphaJump() const;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.5 * 1.3;
};
