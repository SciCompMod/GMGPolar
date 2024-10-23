#include "../include/InputFunctions/ExactSolution/cartesianR2_ShafranovGeometry.h"

CartesianR2_ShafranovGeometry::CartesianR2_ShafranovGeometry(const double& Rmax, const double& elongation_kappa, const double& shift_delta) : 
    Rmax(Rmax),
    elongation_kappa(elongation_kappa),
    shift_delta(shift_delta)
{}

double CartesianR2_ShafranovGeometry::exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (1.0 - (r/Rmax) * (r/Rmax)) * 
        sin(M_PI * (2.0 * elongation_kappa * (r/Rmax) * sin_theta + 2.0 * (r/Rmax) * sin_theta)) * 
        cos(M_PI * ((-2.0) * shift_delta * ((r/Rmax) * (r/Rmax)) - 2.0 * elongation_kappa * (r/Rmax) * cos_theta + 2.0 * (r/Rmax) * cos_theta));
}