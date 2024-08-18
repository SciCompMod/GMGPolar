#include "../include/InputFunctions/BoundaryConditions/cartesianR6_Boundary_ShafranovGeometry.h"

CartesianR6_Boundary_ShafranovGeometry::CartesianR6_Boundary_ShafranovGeometry(const double& Rmax, const double& elongation_kappa, const double& shift_delta) : 
    Rmax(Rmax),
    elongation_kappa(elongation_kappa),
    shift_delta(shift_delta)
{}

double CartesianR6_Boundary_ShafranovGeometry::u_D(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return 0.4096 * pow(((r/Rmax) - 1.0), 6.0) * 
        pow(((r/Rmax) + 1.0), 6.0) * 
        sin(M_PI * (2.0 * elongation_kappa * (r/Rmax) * sin_theta + 2.0 * (r/Rmax) * sin_theta)) * 
        cos(M_PI * ((-2.0) * shift_delta * ((r/Rmax) * (r/Rmax)) - 2.0 * elongation_kappa * (r/Rmax) * cos_theta + 2.0 * (r/Rmax) * cos_theta));
}

double CartesianR6_Boundary_ShafranovGeometry::u_D_Interior(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return 0.4096 * pow(((r/Rmax) - 1.0), 6.0) * 
        pow(((r/Rmax) + 1.0), 6.0) * 
        sin(M_PI * (2.0 * elongation_kappa * (r/Rmax) * sin_theta + 2.0 * (r/Rmax) * sin_theta)) * 
        cos(M_PI * ((-2.0) * shift_delta * ((r/Rmax) * (r/Rmax)) - 2.0 * elongation_kappa * (r/Rmax) * cos_theta + 2.0 * (r/Rmax) * cos_theta));
}
