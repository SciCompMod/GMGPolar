#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class PolarR6_Boundary_CulhamGeometry : public BoundaryConditions
{
public:
    PolarR6_Boundary_CulhamGeometry() = default;
    explicit PolarR6_Boundary_CulhamGeometry(const double& Rmax);
    virtual ~PolarR6_Boundary_CulhamGeometry() = default;

    double u_D(const double& r, const double& theta) const override;
    double u_D_Interior(const double& r, const double& theta) const override;

private:
    const double Rmax = 1.3;
};
