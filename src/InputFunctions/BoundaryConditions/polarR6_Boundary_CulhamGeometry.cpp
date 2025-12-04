#include "../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CulhamGeometry.h"

PolarR6_Boundary_CulhamGeometry::PolarR6_Boundary_CulhamGeometry(double Rmax)
    : Rmax(Rmax)
{
}

double PolarR6_Boundary_CulhamGeometry::u_D(double r, double theta)const
{
    return 0.0;
}

double PolarR6_Boundary_CulhamGeometry::u_D_Interior(double r, double theta)const
{
    return 0.0;
}
