#include "../include/InputFunctions/SourceTerms/polarR6_Poisson_CzarnyGeometry.h"

void PolarR6_Poisson_CzarnyGeometry::initializeGeometry(){
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}

PolarR6_Poisson_CzarnyGeometry::PolarR6_Poisson_CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e) : 
    Rmax(Rmax),
    inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon),
    ellipticity_e(ellipticity_e)
{
    initializeGeometry();
}

double PolarR6_Poisson_CzarnyGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const 
{
    return 0.0;
}
