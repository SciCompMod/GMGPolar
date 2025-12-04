#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class Refined_Boundary_CzarnyGeometry : public BoundaryConditions
{
public:
    explicit Refined_Boundary_CzarnyGeometry();
    explicit Refined_Boundary_CzarnyGeometry(double Rmax, double inverse_aspect_ratio_epsilon, double ellipticity_e);

    virtual ~Refined_Boundary_CzarnyGeometry() = default;

    double u_D(double r, double theta) const override;
    double u_D_Interior(double r, double theta) const override;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
