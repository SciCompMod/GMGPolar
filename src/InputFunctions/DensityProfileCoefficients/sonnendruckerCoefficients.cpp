#include <InputFunctions/DensityProfileCoefficients/sonnendruckerCoefficients.h>
using namespace gmgpolar;

SonnendruckerCoefficients::SonnendruckerCoefficients(double Rmax, double alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

KOKKOS_FUNCTION double SonnendruckerCoefficients::alpha(double r, double theta) const
{
    return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111);
}

KOKKOS_FUNCTION double SonnendruckerCoefficients::beta(double r, double theta) const
{
    return 0.0;
}

KOKKOS_FUNCTION double SonnendruckerCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
