#include "../include/InputFunctions/DensityProfileCoefficients/sonnendruckerGyroCoefficients.h"
using namespace gmgpolar;

SonnendruckerGyroCoefficients::SonnendruckerGyroCoefficients(double _Rmax, double _alpha_jump)
    : Rmax(_Rmax)
    , alpha_jump(_alpha_jump)
{
}

KOKKOS_FUNCTION double SonnendruckerGyroCoefficients::alpha(double r, double theta) const
{
    return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111);
}

KOKKOS_FUNCTION double SonnendruckerGyroCoefficients::beta(double r, double theta) const
{
    return pow((0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)),
               (double)((-1)));
}

KOKKOS_FUNCTION double SonnendruckerGyroCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
