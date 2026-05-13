#include "../include/InputFunctions/DensityProfileCoefficients/zoniCoefficients.h"
using namespace gmgpolar;

ZoniCoefficients::ZoniCoefficients(double Rmax, double alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

KOKKOS_FUNCTION double ZoniCoefficients::alpha(double r, double theta) const
{
    return exp(-tanh(10.0 * (r / Rmax) - 5.0));
}

KOKKOS_FUNCTION double ZoniCoefficients::beta(double r, double theta) const
{
    return 0.0;
}

KOKKOS_FUNCTION double ZoniCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
