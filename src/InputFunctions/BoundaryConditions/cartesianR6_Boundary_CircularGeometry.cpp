#include "../include/InputFunctions/BoundaryConditions/cartesianR6_Boundary_CircularGeometry.h"

CartesianR6_Boundary_CircularGeometry::CartesianR6_Boundary_CircularGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double CartesianR6_Boundary_CircularGeometry::u_D(const double& r, const double& theta)const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return 0.4096 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
           sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta);
}

double CartesianR6_Boundary_CircularGeometry::u_D_Interior(const double& r, const double& theta) const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return 0.4096 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
           sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta);
}
