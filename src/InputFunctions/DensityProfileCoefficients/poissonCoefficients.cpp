#include "../include/InputFunctions/DensityProfileCoefficients/poissonCoefficients.h"

PoissonCoefficients::PoissonCoefficients(double Rmax, double alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

double PoissonCoefficients::alpha(double r, double theta) const
{
    return 1.0;
}

double PoissonCoefficients::beta(double r, double theta) const
{
    return 0.0;
}

double PoissonCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
