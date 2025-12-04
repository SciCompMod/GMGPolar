#include "../include/InputFunctions/DensityProfileCoefficients/zoniCoefficients.h"

ZoniCoefficients::ZoniCoefficients(double Rmax, double alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

double ZoniCoefficients::alpha(double r, double theta) const
{
    return exp(-tanh(10.0 * (r / Rmax) - 5.0));
}

double ZoniCoefficients::beta(double r, double theta) const
{
    return 0.0;
}

double ZoniCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
