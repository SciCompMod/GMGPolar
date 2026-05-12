#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

namespace gmgpolar
{

class ZoniGyroCoefficients
{
public:
    ZoniGyroCoefficients() = default;
    explicit ZoniGyroCoefficients(double Rmax, double alpha);

    double alpha(double r, double theta) const;
    double beta(double r, double theta) const;

    double getAlphaJump() const;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.4837 * 1.3;
};
} // namespace gmgpolar
