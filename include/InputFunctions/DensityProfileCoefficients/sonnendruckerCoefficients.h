#pragma once

#include <cmath>
#include <Kokkos_Core.hpp>

#include "../densityProfileCoefficients.h"

namespace gmgpolar
{

class SonnendruckerCoefficients
{
public:
    SonnendruckerCoefficients() = default;
    explicit SonnendruckerCoefficients(double Rmax, double alpha);

    KOKKOS_FUNCTION double alpha(double r, double theta) const;
    KOKKOS_FUNCTION double beta(double r, double theta) const;

    KOKKOS_FUNCTION double getAlphaJump() const;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.66 * 1.3;
};
} // namespace gmgpolar
