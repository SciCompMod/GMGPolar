#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class PolarR6_Boundary_CulhamGeometry : public BoundaryConditions
{
public:
    PolarR6_Boundary_CulhamGeometry() = default;
    explicit PolarR6_Boundary_CulhamGeometry(double Rmax);
    virtual ~PolarR6_Boundary_CulhamGeometry() = default;

    double u_D(double r, double theta) const override;
    double u_D_Interior(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
