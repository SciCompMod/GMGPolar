#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class Refined_Boundary_CulhamGeometry
{
public:
    explicit Refined_Boundary_CulhamGeometry(double Rmax);

    double u_D(double r, double theta) const;
    double u_D_Interior(double r, double theta) const;

private:
    const double Rmax = 1.3;
};

static_assert(concepts::BoundaryConditions<Refined_Boundary_CulhamGeometry>);
