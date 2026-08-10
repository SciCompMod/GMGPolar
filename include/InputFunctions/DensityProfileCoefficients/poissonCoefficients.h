#pragma once
#include "../../PolarGrid/polargrid.h"

#include <cmath>
#include <Kokkos_Core.hpp>

#include "../densityProfileCoefficients.h"

namespace gmgpolar
{

class PoissonCoefficients
{
public:
    explicit PoissonCoefficients(const PolarGrid& grid, double Rmax, double alpha);

    KOKKOS_FUNCTION double alpha(int i_r, int i_theta) const;
    KOKKOS_FUNCTION double beta(int i_r, int i_theta) const;

    KOKKOS_FUNCTION double getAlphaJump() const;

private:
    PolarGrid grid_;
    const double Rmax       = 1.3;
    const double alpha_jump = 0.5 * 1.3;
};
} // namespace gmgpolar
