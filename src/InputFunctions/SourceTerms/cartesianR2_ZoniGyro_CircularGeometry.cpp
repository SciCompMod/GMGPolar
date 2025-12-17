#include "../include/InputFunctions/SourceTerms/cartesianR2_ZoniGyro_CircularGeometry.h"

CartesianR2_ZoniGyro_CircularGeometry::CartesianR2_ZoniGyro_CircularGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid) , Rmax(Rmax)
{
}

double CartesianR2_ZoniGyro_CircularGeometry::operator()(int i_r, int i_theta) const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return (1.0 - (r / Rmax) * (r / Rmax)) * exp(tanh(10.0 * (r / Rmax) - 5.0)) *
               sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
           ((r / Rmax) * (10.0 * pow(tanh(10.0 * (r / Rmax) - 5.0), 2.0) - 10.0) *
                ((-2.0) * (r / Rmax) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                     cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                 2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) *
                     cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                 2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                     sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta) *
                exp(-tanh(10.0 * (r / Rmax) - 5.0)) +
            (r / Rmax) *
                ((-8.0) * M_PI * (r / Rmax) * sin_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) *
                     cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                 8.0 * M_PI * (r / Rmax) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                     sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta -
                 4.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * pow(sin_theta, 2.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                 8.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                     sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) -
                 4.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                     pow(cos_theta, 2.0) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                 2.0 * sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta)) *
                exp(-tanh(10.0 * (r / Rmax) - 5.0)) +
            ((-2.0) * (r / Rmax) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
             2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) *
                 cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
             2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                 sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta) *
                exp(-tanh(10.0 * (r / Rmax) - 5.0)) +
            ((-4.0) * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * pow(sin_theta, 2.0) *
                 sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
             8.0 * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                 sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) -
             4.0 * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                 sin(2.0 * M_PI * (r / Rmax) * sin_theta) * pow(cos_theta, 2.0) *
                 cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
             2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) *
                 cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
             2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                 sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta) *
                exp(-tanh(10.0 * (r / Rmax) - 5.0))) /
               (r / Rmax);
}