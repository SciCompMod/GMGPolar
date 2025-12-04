#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class PolarR6_Boundary_ShafranovGeometry : public BoundaryConditions
{
public:
    PolarR6_Boundary_ShafranovGeometry() = default;
    explicit PolarR6_Boundary_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);
    virtual ~PolarR6_Boundary_ShafranovGeometry() = default;

    double u_D(double r, double theta) const override;
    double u_D_Interior(double r, double theta) const override;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
