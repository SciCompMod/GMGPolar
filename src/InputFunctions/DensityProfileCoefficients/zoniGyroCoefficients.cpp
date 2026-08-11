#include <InputFunctions/DensityProfileCoefficients/zoniGyroCoefficients.h>
using namespace gmgpolar;

ZoniGyroCoefficients::ZoniGyroCoefficients(const PolarGrid& grid, double Rmax, double alpha_jump)
    : grid_(grid)
    , Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

KOKKOS_FUNCTION double ZoniGyroCoefficients::alpha(int i_r, int i_theta) const
{
    double r = grid_.radius(i_r);
    return exp(-tanh(10.0 * (r / Rmax) - 5.0));
}

KOKKOS_FUNCTION double ZoniGyroCoefficients::beta(int i_r, int i_theta) const
{
    double r = grid_.radius(i_r);
    return exp(tanh(10.0 * (r / Rmax) - 5.0));
}

KOKKOS_FUNCTION double ZoniGyroCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
