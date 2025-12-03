#include "../include/InputFunctions/ExactSolution/polarR6_CircularGeometry.h"

PolarR6_CircularGeometry::PolarR6_CircularGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double PolarR6_CircularGeometry::exact_solution(const double& r, const double& theta)const
{
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
