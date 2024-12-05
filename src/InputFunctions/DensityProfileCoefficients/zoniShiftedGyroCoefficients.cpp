#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedGyroCoefficients.h"

ZoniShiftedGyroCoefficients::ZoniShiftedGyroCoefficients(const double& Rmax, const double& alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

double ZoniShiftedGyroCoefficients::alpha(const double& r) const
{
    return exp(-tanh(20.0 * (r / Rmax) - 14.0));
}

double ZoniShiftedGyroCoefficients::beta(const double& r) const
{
    return exp(tanh(20.0 * (r / Rmax) - 14.0));
}

double ZoniShiftedGyroCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
