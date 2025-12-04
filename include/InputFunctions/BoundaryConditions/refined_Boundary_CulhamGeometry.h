#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class Refined_Boundary_CulhamGeometry : public BoundaryConditions
{
public:
    Refined_Boundary_CulhamGeometry() = default;
    explicit Refined_Boundary_CulhamGeometry(double Rmax);
    virtual ~Refined_Boundary_CulhamGeometry() = default;

    double u_D(double r, double theta) const override;
    double u_D_Interior(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
