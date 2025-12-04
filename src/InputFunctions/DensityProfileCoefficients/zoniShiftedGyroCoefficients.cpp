#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedGyroCoefficients.h"

ZoniShiftedGyroCoefficients::ZoniShiftedGyroCoefficients(double Rmax, double alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

double ZoniShiftedGyroCoefficients::alpha(double r, double theta) const
{
    return exp(-tanh(20.0 * (r / Rmax) - 14.0));
}

double ZoniShiftedGyroCoefficients::beta(double r, double theta) const
{
    return exp(tanh(20.0 * (r / Rmax) - 14.0));
}

double ZoniShiftedGyroCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
