#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class ZoniCoefficients
{
public:
    ZoniCoefficients() = default;
    explicit ZoniCoefficients(double Rmax, double alpha);

    double alpha(double r, double theta) const;
    double beta(double r, double theta) const;

    double getAlphaJump() const;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.4837 * 1.3;
};
