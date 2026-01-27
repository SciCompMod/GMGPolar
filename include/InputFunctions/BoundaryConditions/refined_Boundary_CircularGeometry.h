#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class Refined_Boundary_CircularGeometry
{
public:
    Refined_Boundary_CircularGeometry() = default;
    explicit Refined_Boundary_CircularGeometry(double Rmax);

    double u_D(double r, double theta) const;
    double u_D_Interior(double r, double theta) const;

private:
    const double Rmax = 1.3;
};

static_assert(concepts::BoundaryConditions<Refined_Boundary_CircularGeometry>);
