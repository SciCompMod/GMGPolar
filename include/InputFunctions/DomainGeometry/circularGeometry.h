#pragma once

#include <cmath>
#include <Kokkos_Core.hpp>

#include "../domainGeometry.h"

namespace gmgpolar
{

class CircularGeometry
{
public:
    CircularGeometry() = default;
    explicit CircularGeometry(double Rmax);

    KOKKOS_INLINE_FUNCTION double Fx(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double Fy(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFx_dr(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFy_dr(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFx_dtheta(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFy_dtheta(double r, double theta) const;

private:
    const double Rmax = 1.3;
};

#include "circularGeometry.inl"
} // namespace gmgpolar
