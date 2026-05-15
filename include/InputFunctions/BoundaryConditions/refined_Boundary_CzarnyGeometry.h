#pragma once

#include <cmath>

#include "../boundaryConditions.h"

namespace gmgpolar
{

class Refined_Boundary_CzarnyGeometry
{
public:
    explicit Refined_Boundary_CzarnyGeometry();
    explicit Refined_Boundary_CzarnyGeometry(double Rmax, double inverse_aspect_ratio_epsilon, double ellipticity_e);
KOKKOS_DEFAULTED_FUNCTION Refined_Boundary_CzarnyGeometry(const Refined_Boundary_CzarnyGeometry&) = default;

    KOKKOS_FUNCTION double u_D(double r, double theta) const;
    KOKKOS_FUNCTION double u_D_Interior(double r, double theta) const;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};

static_assert(concepts::BoundaryConditions<Refined_Boundary_CzarnyGeometry>);
} // namespace gmgpolar
