#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CulhamGeometry.h"

PolarR6_ZoniShiftedGyro_CulhamGeometry::PolarR6_ZoniShiftedGyro_CulhamGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double PolarR6_ZoniShiftedGyro_CulhamGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta,
                                                     const double& cos_theta) const
{
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}