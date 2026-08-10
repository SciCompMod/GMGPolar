#include <InputFunctions/DensityProfileCoefficients/zoniShiftedGyroCoefficients.h>
using namespace gmgpolar;

ZoniShiftedGyroCoefficients::ZoniShiftedGyroCoefficients(const PolarGrid& grid, double Rmax, double alpha_jump)
    : grid_(grid)
    , Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

KOKKOS_FUNCTION double ZoniShiftedGyroCoefficients::alpha(int i_r, int i_theta) const
{
    double r = grid_.radius(i_r);
    return exp(-tanh(20.0 * (r / Rmax) - 14.0));
}

KOKKOS_FUNCTION double ZoniShiftedGyroCoefficients::beta(int i_r, int i_theta) const
{
    double r = grid_.radius(i_r);
    return exp(tanh(20.0 * (r / Rmax) - 14.0));
}

KOKKOS_FUNCTION double ZoniShiftedGyroCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
