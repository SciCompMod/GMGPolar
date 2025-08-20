#include "../include/InputFunctions/DensityProfileCoefficients/zoniGyroCoefficients.h"

ZoniGyroCoefficients::ZoniGyroCoefficients(const double& Rmax, const double& alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

double ZoniGyroCoefficients::alpha(const double& r, const double& theta) const
{
    return exp(-tanh(10.0 * (r / Rmax) - 5.0));
}

double ZoniGyroCoefficients::beta(const double& r, const double& theta) const
{
    return exp(tanh(10.0 * (r / Rmax) - 5.0));
}

double ZoniGyroCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
