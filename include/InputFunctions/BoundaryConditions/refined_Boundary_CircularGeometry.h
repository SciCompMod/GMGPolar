#pragma once

#include <cmath>
#include <Kokkos_Core.hpp>

#include "../boundaryConditions.h"

namespace gmgpolar
{

class Refined_Boundary_CircularGeometry
{
public:
    explicit Refined_Boundary_CircularGeometry(double Rmax);
    KOKKOS_DEFAULTED_FUNCTION Refined_Boundary_CircularGeometry(const Refined_Boundary_CircularGeometry&) = default;

    KOKKOS_FUNCTION double u_D(double r, double theta) const;
    KOKKOS_FUNCTION double u_D_Interior(double r, double theta) const;

private:
    const double Rmax = 1.3;
};

static_assert(concepts::BoundaryConditions<Refined_Boundary_CircularGeometry>);
} // namespace gmgpolar
