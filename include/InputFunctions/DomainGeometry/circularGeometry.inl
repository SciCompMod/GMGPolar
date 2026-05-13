#pragma once

#include "circularGeometry.h"

// In earlier versions denoted by 'x'
KOKKOS_INLINE_FUNCTION double CircularGeometry::Fx(double r, double theta) const
{
    return (r / Rmax) * std::cos(theta);
}

// In earlier versions denoted by 'y'
KOKKOS_INLINE_FUNCTION double CircularGeometry::Fy(double r, double theta) const
{
    return (r / Rmax) * std::sin(theta);
}

// In earlier versions denoted by 'Jrr'
KOKKOS_INLINE_FUNCTION double CircularGeometry::dFx_dr(double r, double theta) const
{
    return (std::cos(theta)) / Rmax;
}

// In earlier versions denoted by 'Jtr'
KOKKOS_INLINE_FUNCTION double CircularGeometry::dFy_dr(double r, double theta) const
{
    return (std::sin(theta)) / Rmax;
}

// In earlier versions denoted by 'Jrt'
KOKKOS_INLINE_FUNCTION double CircularGeometry::dFx_dt(double r, double theta) const
{
    return (-(r / Rmax)) * std::sin(theta);
}

// In earlier versions denoted by 'Jtt'
KOKKOS_INLINE_FUNCTION double CircularGeometry::dFy_dt(double r, double theta) const
{
    return (r / Rmax) * std::cos(theta);
}
