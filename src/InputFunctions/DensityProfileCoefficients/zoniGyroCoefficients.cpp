#include "../include/InputFunctions/DensityProfileCoefficients/zoniGyroCoefficients.h"

ZoniGyroCoefficients::ZoniGyroCoefficients(double Rmax, double alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

double ZoniGyroCoefficients::alpha(double r, double theta) const
{
    return exp(-tanh(10.0 * (r / Rmax) - 5.0));
}

double ZoniGyroCoefficients::beta(double r, double theta) const
{
    return exp(tanh(10.0 * (r / Rmax) - 5.0));
}

double ZoniGyroCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
