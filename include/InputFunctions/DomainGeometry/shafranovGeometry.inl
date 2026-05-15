#pragma once

#include "shafranovGeometry.h"

// In earlier versions denoted by 'x'
KOKKOS_INLINE_FUNCTION double ShafranovGeometry::Fx(double r, double theta) const
{
    double cos_theta = std::cos(theta);
    return (1.0 - elongation_kappa) * (r / Rmax) * cos_theta - shift_delta * (r / Rmax) * (r / Rmax);
}

// In earlier versions denoted by 'y'
KOKKOS_INLINE_FUNCTION double ShafranovGeometry::Fy(double r, double theta) const
{
    double sin_theta = std::sin(theta);
    return (1.0 + elongation_kappa) * (r / Rmax) * sin_theta;
}

// In earlier versions denoted by 'Jrr'
KOKKOS_INLINE_FUNCTION double ShafranovGeometry::dFx_dr(double r, double theta) const
{
    double cos_theta = std::cos(theta);
    return ((Rmax - elongation_kappa * Rmax) * cos_theta - 2.0 * shift_delta * r) / (Rmax * Rmax);
}

// In earlier versions denoted by 'Jtr'
KOKKOS_INLINE_FUNCTION double ShafranovGeometry::dFy_dr(double r, double theta) const
{
    double sin_theta = std::sin(theta);
    return (elongation_kappa + 1.0) * sin_theta / Rmax;
}

// In earlier versions denoted by 'Jrt'
KOKKOS_INLINE_FUNCTION double ShafranovGeometry::dFx_dt(double r, double theta) const
{
    double sin_theta = std::sin(theta);
    return ((elongation_kappa - 1.0) * r * sin_theta) / Rmax;
}

// In earlier versions denoted by 'Jtt'
KOKKOS_INLINE_FUNCTION double ShafranovGeometry::dFy_dt(double r, double theta) const
{
    double cos_theta = std::cos(theta);
    return ((elongation_kappa + 1.0) * r * cos_theta) / Rmax;
}
