#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h"
#include <iostream>
using namespace gmgpolar;
ZoniShiftedCoefficients::ZoniShiftedCoefficients(double Rmax, double alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

KOKKOS_FUNCTION double ZoniShiftedCoefficients::alpha(double r, double theta) const
{
    return exp(-tanh(20.0 * (r / Rmax) - 14.0));
}

KOKKOS_FUNCTION double ZoniShiftedCoefficients::beta(double r, double theta) const
{
    return 0.0;
}

KOKKOS_FUNCTION double ZoniShiftedCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
