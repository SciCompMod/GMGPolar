#include "../include/InputFunctions/BoundaryConditions/polarR6_Boundary_ShafranovGeometry.h"

PolarR6_Boundary_ShafranovGeometry::PolarR6_Boundary_ShafranovGeometry(double Rmax, double elongation_kappa,
                                                                       double shift_delta)
    : Rmax(Rmax)
    , elongation_kappa(elongation_kappa)
    , shift_delta(shift_delta)
{
}

double PolarR6_Boundary_ShafranovGeometry::u_D(double r, double theta) const
{
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}

double PolarR6_Boundary_ShafranovGeometry::u_D_Interior(double r, double theta) const
{
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
