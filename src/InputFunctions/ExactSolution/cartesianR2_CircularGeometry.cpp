#include "../include/InputFunctions/ExactSolution/cartesianR2_CircularGeometry.h"

CartesianR2_CircularGeometry::CartesianR2_CircularGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double CartesianR2_CircularGeometry::exact_solution(const double& r, const double& theta) const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
           cos(2.0 * M_PI * (r / Rmax) * cos_theta);
}