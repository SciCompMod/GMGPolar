#include "../include/InputFunctions/SourceTerms/cartesianR2_Poisson_ShafranovGeometry.h"

CartesianR2_Poisson_ShafranovGeometry::CartesianR2_Poisson_ShafranovGeometry(const double& Rmax, const double& elongation_kappa, const double& shift_delta) : 
    Rmax(Rmax),
    elongation_kappa(elongation_kappa),
    shift_delta(shift_delta)
{}

double CartesianR2_Poisson_ShafranovGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const 
{
    return 0.0;
}