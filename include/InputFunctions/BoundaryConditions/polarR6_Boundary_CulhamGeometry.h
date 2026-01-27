#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class PolarR6_Boundary_CulhamGeometry
{
public:
    explicit PolarR6_Boundary_CulhamGeometry(double Rmax);

    double u_D(double r, double theta) const;
    double u_D_Interior(double r, double theta) const;

private:
    const double Rmax = 1.3;
};

static_assert(concepts::BoundaryConditions<PolarR6_Boundary_CulhamGeometry>);
