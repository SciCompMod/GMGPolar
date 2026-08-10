#include <InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h>
#include <iostream>
using namespace gmgpolar;
ZoniShiftedCoefficients::ZoniShiftedCoefficients(const PolarGrid& grid, double Rmax, double alpha_jump)
    : grid_(grid)
    , Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

KOKKOS_FUNCTION double ZoniShiftedCoefficients::alpha(int i_r, int i_theta) const
{
    double r = grid_.radius(i_r);
    return exp(-tanh(20.0 * (r / Rmax) - 14.0));
}

KOKKOS_FUNCTION double ZoniShiftedCoefficients::beta(int i_r, int i_theta) const
{
    return 0.0;
}

KOKKOS_FUNCTION double ZoniShiftedCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
