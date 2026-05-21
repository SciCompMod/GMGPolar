#include <InputFunctions/BoundaryConditions/polarR6_Boundary_CircularGeometry.h>
using namespace gmgpolar;

PolarR6_Boundary_CircularGeometry::PolarR6_Boundary_CircularGeometry(double Rmax)
    : Rmax(Rmax)
{
}

KOKKOS_FUNCTION double PolarR6_Boundary_CircularGeometry::u_D(double r, double theta) const
{
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}

KOKKOS_FUNCTION double PolarR6_Boundary_CircularGeometry::u_D_Interior(double r, double theta) const
{
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
