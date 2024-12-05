#include "../include/InputFunctions/ExactSolution/polarR6_CulhamGeometry.h"

PolarR6_CulhamGeometry::PolarR6_CulhamGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double PolarR6_CulhamGeometry::exact_solution(const double& r, const double& theta, const double& sin_theta,
                                              const double& cos_theta) const
{
    return 0.0;
}
