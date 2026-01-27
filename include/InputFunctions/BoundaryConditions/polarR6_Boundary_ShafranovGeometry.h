#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class PolarR6_Boundary_ShafranovGeometry
{
public:
    explicit PolarR6_Boundary_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);

    double u_D(double r, double theta) const;
    double u_D_Interior(double r, double theta) const;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};

static_assert(concepts::BoundaryConditions<PolarR6_Boundary_ShafranovGeometry>);
