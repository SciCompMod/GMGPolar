#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class CartesianR6_Boundary_CircularGeometry
{
public:
    CartesianR6_Boundary_CircularGeometry() = default;
    explicit CartesianR6_Boundary_CircularGeometry(double Rmax);
    virtual ~CartesianR6_Boundary_CircularGeometry() = default;

    double u_D(double r, double theta) const ;
    double u_D_Interior(double r, double theta) const ;

private:
    const double Rmax = 1.3;
};

static_assert(concepts::BoundaryConditions<CartesianR6_Boundary_CircularGeometry>);
