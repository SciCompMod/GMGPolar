#include "../include/InputFunctions/DensityProfileCoefficients/zoniCoefficients.h"

ZoniCoefficients::ZoniCoefficients(const double& Rmax, const double& alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

double ZoniCoefficients::alpha(const double& r, const double& theta) const
{
    return exp(-tanh(10.0 * (r / Rmax) - 5.0));
}

double ZoniCoefficients::beta(const double& r, const double& theta) const
{
    return 0.0;
}

double ZoniCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
