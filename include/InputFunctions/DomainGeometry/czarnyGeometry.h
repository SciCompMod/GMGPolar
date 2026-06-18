#pragma once

#include <cmath>
#include <Kokkos_Core.hpp>

#include "../domainGeometry.h"

/* Triangular Shape and Ellipticity */

namespace gmgpolar
{

class CzarnyGeometry
{
public:
    explicit CzarnyGeometry();
    explicit CzarnyGeometry(double Rmax, double inverse_aspect_ratio_epsilon, double ellipticity_e);

    KOKKOS_INLINE_FUNCTION double Fx(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double Fy(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFx_dr(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFy_dr(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFx_dtheta(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFy_dtheta(double r, double theta) const;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};

#include "czarnyGeometry.inl"
} // namespace gmgpolar
