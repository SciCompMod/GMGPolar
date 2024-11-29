#pragma once

#include <cmath>

#include "../boundaryConditions.h"

class Refined_Boundary_CircularGeometry : public BoundaryConditions { 
public:
    Refined_Boundary_CircularGeometry() = default;
    explicit Refined_Boundary_CircularGeometry(const double& Rmax);
    virtual ~Refined_Boundary_CircularGeometry() = default;

    double u_D(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double u_D_Interior(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    
private:
    const double Rmax = 1.3;
};
