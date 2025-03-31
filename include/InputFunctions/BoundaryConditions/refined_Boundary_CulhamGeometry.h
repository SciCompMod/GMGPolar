#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class Refined_Boundary_CulhamGeometry : public BoundaryConditions
{
public:
    Refined_Boundary_CulhamGeometry() = default;
    explicit Refined_Boundary_CulhamGeometry(const double& Rmax);
    virtual ~Refined_Boundary_CulhamGeometry() = default;

    double u_D(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double u_D_Interior(const double& r, const double& theta, const double& sin_theta,
                        const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};
