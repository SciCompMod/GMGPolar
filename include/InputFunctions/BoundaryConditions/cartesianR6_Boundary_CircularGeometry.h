#pragma once

#include <cmath>

#include "../boundaryConditions.h"

namespace gmgpolar
{

class CartesianR6_Boundary_CircularGeometry
{
public:
    explicit CartesianR6_Boundary_CircularGeometry(double Rmax);
KOKKOS_DEFAULTED_FUNCTION CartesianR6_Boundary_CircularGeometry(const CartesianR6_Boundary_CircularGeometry&) = default;

    KOKKOS_FUNCTION double u_D(double r, double theta) const;
    KOKKOS_FUNCTION double u_D_Interior(double r, double theta) const;

private:
    const double Rmax = 1.3;
};

static_assert(concepts::BoundaryConditions<CartesianR6_Boundary_CircularGeometry>);
} // namespace gmgpolar
