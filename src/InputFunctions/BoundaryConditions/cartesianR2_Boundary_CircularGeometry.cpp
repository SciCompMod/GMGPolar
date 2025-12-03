#include "../include/InputFunctions/BoundaryConditions/cartesianR2_Boundary_CircularGeometry.h"

CartesianR2_Boundary_CircularGeometry::CartesianR2_Boundary_CircularGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double CartesianR2_Boundary_CircularGeometry::u_D(const double& r, const double& theta)const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
           cos(2.0 * M_PI * (r / Rmax) * cos_theta);
}

double CartesianR2_Boundary_CircularGeometry::u_D_Interior(const double& r, const double& theta,
                                                           const double& sin_theta, const double& cos_theta) const
{
    return (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
           cos(2.0 * M_PI * (r / Rmax) * cos_theta);
}