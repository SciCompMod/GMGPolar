#include "../include/InputFunctions/DensityProfileCoefficients/poissonCoefficients.h"

PoissonCoefficients::PoissonCoefficients(const double& Rmax, const double& alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

double PoissonCoefficients::alpha(const double& r) const
{
    return 1.0;
}

double PoissonCoefficients::beta(const double& r) const
{
    return 0.0;
}

double PoissonCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
