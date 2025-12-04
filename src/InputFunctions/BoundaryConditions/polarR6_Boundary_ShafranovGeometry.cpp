#include "../include/InputFunctions/BoundaryConditions/polarR6_Boundary_ShafranovGeometry.h"

PolarR6_Boundary_ShafranovGeometry::PolarR6_Boundary_ShafranovGeometry(const double& Rmax,
                                                                       const double& elongation_kappa,
                                                                       const double& shift_delta)
    : Rmax(Rmax)
    , elongation_kappa(elongation_kappa)
    , shift_delta(shift_delta)
{
}

double PolarR6_Boundary_ShafranovGeometry::u_D(const double& r, const double& theta) const
{
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}

double PolarR6_Boundary_ShafranovGeometry::u_D_Interior(const double& r, const double& theta) const
{
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
