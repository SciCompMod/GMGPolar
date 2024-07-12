#include "../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CircularGeometry.h"

PolarR6_Boundary_CircularGeometry::PolarR6_Boundary_CircularGeometry(const double& Rmax) : 
    Rmax(Rmax) 
{}

double PolarR6_Boundary_CircularGeometry::u_D(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}

double PolarR6_Boundary_CircularGeometry::u_D_Interior(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
