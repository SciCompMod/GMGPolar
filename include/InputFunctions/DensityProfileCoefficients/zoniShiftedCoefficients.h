#pragma once

#include <cmath>
#include <Kokkos_Core.hpp>

#include "../densityProfileCoefficients.h"

namespace gmgpolar
{

class ZoniShiftedCoefficients
{
public:
    ZoniShiftedCoefficients() = default;
    explicit ZoniShiftedCoefficients(double Rmax, double alpha);

    KOKKOS_FUNCTION double alpha(double r, double theta) const;
    KOKKOS_FUNCTION double beta(double r, double theta) const;

    KOKKOS_FUNCTION double getAlphaJump() const;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.678 * 1.3;
};
} // namespace gmgpolar
