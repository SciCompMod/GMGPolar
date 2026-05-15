#pragma once

#include <cmath>
#include <Kokkos_Core.hpp>

#include "../boundaryConditions.h"

namespace gmgpolar
{

class PolarR6_Boundary_CulhamGeometry
{
public:
    explicit PolarR6_Boundary_CulhamGeometry(double Rmax);
KOKKOS_DEFAULTED_FUNCTION PolarR6_Boundary_CulhamGeometry(const PolarR6_Boundary_CulhamGeometry&) = default;

    KOKKOS_FUNCTION double u_D(double r, double theta) const;
    KOKKOS_FUNCTION double u_D_Interior(double r, double theta) const;

private:
    const double Rmax = 1.3;
};

static_assert(concepts::BoundaryConditions<PolarR6_Boundary_CulhamGeometry>);
} // namespace gmgpolar
