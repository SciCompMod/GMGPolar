#include "../include/InputFunctions/BoundaryConditions/refined_Boundary_CulhamGeometry.h"

Refined_Boundary_CulhamGeometry::Refined_Boundary_CulhamGeometry(double Rmax)
    : Rmax(Rmax)
{
}

double Refined_Boundary_CulhamGeometry::u_D(double r, double theta) const
{
    return 0.0;
}

double Refined_Boundary_CulhamGeometry::u_D_Interior(double r, double theta) const
{
    return 0.0;
}
