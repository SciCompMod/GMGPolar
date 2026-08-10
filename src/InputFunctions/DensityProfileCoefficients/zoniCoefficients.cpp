#include <InputFunctions/DensityProfileCoefficients/zoniCoefficients.h>
using namespace gmgpolar;

ZoniCoefficients::ZoniCoefficients(const PolarGrid& grid, double Rmax, double alpha_jump)
    : grid_(grid)
    , Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

KOKKOS_FUNCTION double ZoniCoefficients::alpha(int i_r, int i_theta) const
{
    double r = grid_.radius(i_r);
    return exp(-tanh(10.0 * (r / Rmax) - 5.0));
}

KOKKOS_FUNCTION double ZoniCoefficients::beta(int i_r, int i_theta) const
{
    return 0.0;
}

KOKKOS_FUNCTION double ZoniCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
