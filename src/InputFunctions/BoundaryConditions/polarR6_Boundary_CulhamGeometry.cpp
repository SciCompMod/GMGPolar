#include "../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CulhamGeometry.h"

PolarR6_Boundary_CulhamGeometry::PolarR6_Boundary_CulhamGeometry(const double& Rmax) : 
    Rmax(Rmax) 
{}

double PolarR6_Boundary_CulhamGeometry::u_D(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return 0.0;
}

double PolarR6_Boundary_CulhamGeometry::u_D_Interior(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return 0.0;
}