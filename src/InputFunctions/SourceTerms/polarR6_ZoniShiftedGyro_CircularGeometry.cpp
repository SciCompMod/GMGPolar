#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CircularGeometry.h"

PolarR6_ZoniShiftedGyro_CircularGeometry::PolarR6_ZoniShiftedGyro_CircularGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double PolarR6_ZoniShiftedGyro_CircularGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta,
                                                       const double& cos_theta) const
{
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * exp(tanh(20.0 * (r / Rmax) - 14.0)) *
               cos(11.0 * theta) -
           pow((r / Rmax), 4.0) *
               ((r / Rmax) *
                    (12.288 * (r / Rmax) * pow(((r / Rmax) - 1.0), 4.0) * cos(11.0 * theta) +
                     17.2032 * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta)) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) +
                (r / Rmax) *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                    (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r / Rmax) - 14.0)) -
                49.5616 * pow(((r / Rmax) - 1.0), 6.0) * exp(-tanh(20.0 * (r / Rmax) - 14.0)) * cos(11.0 * theta) +
                6.0 *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)));
}