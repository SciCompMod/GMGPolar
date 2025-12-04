#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class PolarR6_Boundary_CzarnyGeometry : public BoundaryConditions
{
public:
    explicit PolarR6_Boundary_CzarnyGeometry();
    explicit PolarR6_Boundary_CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon,
                                             const double& ellipticity_e);

    virtual ~PolarR6_Boundary_CzarnyGeometry() = default;

    double u_D(const double& r, const double& theta) const override;
    double u_D_Interior(const double& r, const double& theta) const override;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
