#include <InputFunctions/BoundaryConditions/refined_Boundary_CulhamGeometry.h>
using namespace gmgpolar;

Refined_Boundary_CulhamGeometry::Refined_Boundary_CulhamGeometry(double Rmax)
    : Rmax(Rmax)
{
}

KOKKOS_FUNCTION double Refined_Boundary_CulhamGeometry::u_D(double r, double theta) const
{
    return 0.0;
}

KOKKOS_FUNCTION double Refined_Boundary_CulhamGeometry::u_D_Interior(double r, double theta) const
{
    return 0.0;
}
