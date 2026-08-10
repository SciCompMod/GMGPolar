#include <InputFunctions/DensityProfileCoefficients/sonnendruckerCoefficients.h>
using namespace gmgpolar;

SonnendruckerCoefficients::SonnendruckerCoefficients(const PolarGrid& grid, double Rmax, double alpha_jump)
    : grid_(grid)
    , Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

KOKKOS_FUNCTION double SonnendruckerCoefficients::alpha(int i_r, int i_theta) const
{
    double r = grid_.radius(i_r);
    return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111);
}

KOKKOS_FUNCTION double SonnendruckerCoefficients::beta(int i_r, int i_theta) const
{
    return 0.0;
}

KOKKOS_FUNCTION double SonnendruckerCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
