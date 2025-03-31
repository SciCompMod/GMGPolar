#include "../include/InputFunctions/BoundaryConditions/refined_Boundary_CulhamGeometry.h"

Refined_Boundary_CulhamGeometry::Refined_Boundary_CulhamGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double Refined_Boundary_CulhamGeometry::u_D(const double& r, const double& theta, const double& sin_theta,
                                            const double& cos_theta) const
{
    return 0.0;
}

double Refined_Boundary_CulhamGeometry::u_D_Interior(const double& r, const double& theta, const double& sin_theta,
                                                     const double& cos_theta) const
{
    return 0.0;
}
