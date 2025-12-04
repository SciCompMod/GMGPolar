#include "../include/InputFunctions/ExactSolution/cartesianR6_CzarnyGeometry.h"

void CartesianR6_CzarnyGeometry::initializeGeometry()
{
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}

CartesianR6_CzarnyGeometry::CartesianR6_CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon,
                                                       const double& ellipticity_e)
    : Rmax(Rmax)
    , inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon)
    , ellipticity_e(ellipticity_e)
{
    initializeGeometry();
}

double CartesianR6_CzarnyGeometry::exact_solution(const double& r, const double& theta) const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double temp =
        sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r / Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0);
    return 0.4096 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
           sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / (2.0 - temp)) *
           cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon);
}
