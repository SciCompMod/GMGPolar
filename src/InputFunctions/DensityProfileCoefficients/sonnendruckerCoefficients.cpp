#include "../include/InputFunctions/DensityProfileCoefficients/sonnendruckerCoefficients.h"

SonnendruckerCoefficients::SonnendruckerCoefficients(const double& Rmax, const double& alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

double SonnendruckerCoefficients::alpha(const double& r) const
{
    return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111);
}

double SonnendruckerCoefficients::beta(const double& r) const
{
    return 0.0;
}

double SonnendruckerCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
