#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedGyroCoefficients.h"
using namespace gmgpolar;

ZoniShiftedGyroCoefficients::ZoniShiftedGyroCoefficients(double Rmax, double alpha_jump)
    : Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

KOKKOS_FUNCTION double ZoniShiftedGyroCoefficients::alpha(double r, double theta) const
{
    return exp(-tanh(20.0 * (r / Rmax) - 14.0));
}

KOKKOS_FUNCTION double ZoniShiftedGyroCoefficients::beta(double r, double theta) const
{
    return exp(tanh(20.0 * (r / Rmax) - 14.0));
}

KOKKOS_FUNCTION double ZoniShiftedGyroCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
