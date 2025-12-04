#include "../include/InputFunctions/ExactSolution/cartesianR6_CircularGeometry.h"

CartesianR6_CircularGeometry::CartesianR6_CircularGeometry(double Rmax)
    : Rmax(Rmax)
{
}

double CartesianR6_CircularGeometry::exact_solution(double r, double theta)const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return 0.4096 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
           sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta);
}