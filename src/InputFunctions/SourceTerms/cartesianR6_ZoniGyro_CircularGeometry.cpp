#include "../include/InputFunctions/SourceTerms/cartesianR6_ZoniGyro_CircularGeometry.h"

CartesianR6_ZoniGyro_CircularGeometry::CartesianR6_ZoniGyro_CircularGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double CartesianR6_ZoniGyro_CircularGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta,
                                                    const double& cos_theta) const
{
    return 0.4096 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * exp(tanh(10.0 * (r / Rmax) - 5.0)) *
               sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
           ((r / Rmax) * (10.0 * pow(tanh(10.0 * (r / Rmax) - 5.0), 2.0) - 10.0) *
                (0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                     cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                 0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 5.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 6.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta)) *
                exp(-tanh(10.0 * (r / Rmax) - 5.0)) +
            (r / Rmax) *
                ((-1.6384) * (M_PI * M_PI) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                     pow(sin_theta, 2.0) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                     cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                 3.2768 * (M_PI * M_PI) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                     sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) -
                 1.6384 * (M_PI * M_PI) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * pow(cos_theta, 2.0) *
                     cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                 9.8304 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 5.0) * sin_theta *
                     cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                 9.8304 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 5.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
                 12.288 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 4.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                 9.8304 * M_PI * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                     cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                 9.8304 * M_PI * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 6.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
                 29.4912 * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 5.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                 12.288 * pow(((r / Rmax) - 1.0), 4.0) * pow(((r / Rmax) + 1.0), 6.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta)) *
                exp(-tanh(10.0 * (r / Rmax) - 5.0)) +
            (0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                 cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
             0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                 sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
             2.4576 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 5.0) *
                 sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
             2.4576 * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 6.0) *
                 sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta)) *
                exp(-tanh(10.0 * (r / Rmax) - 5.0)) +
            ((-1.6384) * (M_PI * M_PI) * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                 pow(sin_theta, 2.0) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                 cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
             3.2768 * (M_PI * M_PI) * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                 sin_theta * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta *
                 cos(2.0 * M_PI * (r / Rmax) * sin_theta) -
             1.6384 * (M_PI * M_PI) * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                 sin(2.0 * M_PI * (r / Rmax) * sin_theta) * pow(cos_theta, 2.0) *
                 cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
             0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                 cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
             0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                 sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta) *
                exp(-tanh(10.0 * (r / Rmax) - 5.0))) /
               (r / Rmax);
}