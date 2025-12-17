#include "../include/InputFunctions/SourceTerms/cartesianR2_ZoniShifted_CircularGeometry.h"

CartesianR2_ZoniShifted_CircularGeometry::CartesianR2_ZoniShifted_CircularGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid) , Rmax(Rmax)
{
}

double CartesianR2_ZoniShifted_CircularGeometry::operator()(int i_r, int i_theta) const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return (-((r / Rmax) * (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) *
                  ((-2.0) * (r / Rmax) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                       cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) *
                       cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                       sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) +
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
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) +
              ((-2.0) * (r / Rmax) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                   cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
               2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) *
                   cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
               2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                   sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) +
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
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)))) /
           (r / Rmax);
}