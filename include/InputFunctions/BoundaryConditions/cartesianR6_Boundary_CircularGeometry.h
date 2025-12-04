#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class CartesianR6_Boundary_CircularGeometry : public BoundaryConditions
{
public:
    CartesianR6_Boundary_CircularGeometry() = default;
    explicit CartesianR6_Boundary_CircularGeometry(double Rmax);
    virtual ~CartesianR6_Boundary_CircularGeometry() = default;

    double u_D(double r, double theta) const override;
    double u_D_Interior(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
