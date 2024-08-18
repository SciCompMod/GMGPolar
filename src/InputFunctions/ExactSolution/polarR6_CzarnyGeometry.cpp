#include "../include/InputFunctions/ExactSolution/polarR6_CzarnyGeometry.h"

void PolarR6_CzarnyGeometry::initializeGeometry(){
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}

PolarR6_CzarnyGeometry::PolarR6_CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e) : 
    Rmax(Rmax),
    inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon),
    ellipticity_e(ellipticity_e)
{
    initializeGeometry();
}

double PolarR6_CzarnyGeometry::exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}

