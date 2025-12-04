#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class CartesianR2_Boundary_CircularGeometry : public BoundaryConditions
{
public:
    CartesianR2_Boundary_CircularGeometry() = default;
    explicit CartesianR2_Boundary_CircularGeometry(const double& Rmax);
    virtual ~CartesianR2_Boundary_CircularGeometry() = default;

    double u_D(const double& r, const double& theta) const override;
    double u_D_Interior(const double& r, const double& theta) const override;

private:
    const double Rmax = 1.3;
};
