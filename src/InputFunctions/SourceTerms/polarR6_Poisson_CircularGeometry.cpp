#include "../include/InputFunctions/SourceTerms/polarR6_Poisson_CircularGeometry.h"

PolarR6_Poisson_CircularGeometry::PolarR6_Poisson_CircularGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double PolarR6_Poisson_CircularGeometry::rhs_f(const double& r, const double& theta)const
{
    return (-pow((r / Rmax), 4.0)) * (14.7456 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                                      1.0 * (r / Rmax) *
                                          (12.288 * (r / Rmax) * pow(((r / Rmax) - 1.0), 4.0) * cos(11.0 * theta) +
                                           17.2032 * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta)) -
                                      34.816 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta));
}
