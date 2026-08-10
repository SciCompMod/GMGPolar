#include <InputFunctions/DensityProfileCoefficients/poissonCoefficients.h>
using namespace gmgpolar;

PoissonCoefficients::PoissonCoefficients(const PolarGrid& grid, double Rmax, double alpha_jump)
    : grid_(grid)
    , Rmax(Rmax)
    , alpha_jump(alpha_jump)
{
}

KOKKOS_FUNCTION double PoissonCoefficients::alpha(double r, double theta) const
{
    return 1.0;
}

KOKKOS_FUNCTION double PoissonCoefficients::beta(double r, double theta) const
{
    return 0.0;
}

KOKKOS_FUNCTION double PoissonCoefficients::getAlphaJump() const
{
    return alpha_jump;
}
