#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class CartesianR6_Boundary_ShafranovGeometry : public BoundaryConditions
{
public:
    CartesianR6_Boundary_ShafranovGeometry() = default;
    explicit CartesianR6_Boundary_ShafranovGeometry(const double& Rmax, const double& elongation_kappa,
                                                    const double& shift_delta);
    virtual ~CartesianR6_Boundary_ShafranovGeometry() = default;

    double u_D(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double u_D_Interior(const double& r, const double& theta, const double& sin_theta,
                        const double& cos_theta) const override;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
