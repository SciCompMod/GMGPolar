#pragma once

#include <cmath>
#include <Kokkos_Core.hpp>

#include "../boundaryConditions.h"

namespace gmgpolar
{

class Refined_Boundary_ShafranovGeometry
{
public:
    explicit Refined_Boundary_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);
    KOKKOS_DEFAULTED_FUNCTION Refined_Boundary_ShafranovGeometry(const Refined_Boundary_ShafranovGeometry&) = default;

    KOKKOS_FUNCTION double u_D(double r, double theta) const;
    KOKKOS_FUNCTION double u_D_Interior(double r, double theta) const;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};

static_assert(concepts::BoundaryConditions<Refined_Boundary_ShafranovGeometry>);
} // namespace gmgpolar
