#include "../include/InputFunctions/SourceTerms/cartesianR6_ZoniShiftedGyro_CircularGeometry.h"

CartesianR6_ZoniShiftedGyro_CircularGeometry::CartesianR6_ZoniShiftedGyro_CircularGeometry(PolarGrid const& grid,
                                                                                           double Rmax)
    : grid_(grid)
    , Rmax(Rmax)
{
}

double CartesianR6_ZoniShiftedGyro_CircularGeometry::operator()(std::size_t i_r, std::size_t i_theta) const
{
    double r         = grid_.radius(i_r);
    double theta     = grid_.theta(i_theta);
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return 0.4096 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * exp(tanh(20.0 * (r / Rmax) - 14.0)) *
               sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
           ((r / Rmax) * (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) *
                (0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                     cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                 0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 5.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 6.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) +
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
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) +
            (0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                 cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
             0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                 sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
             2.4576 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 5.0) *
                 sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
             2.4576 * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 6.0) *
                 sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) +
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
                exp(-tanh(20.0 * (r / Rmax) - 14.0))) /
               (r / Rmax);
}