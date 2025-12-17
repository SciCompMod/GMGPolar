#include "../include/InputFunctions/SourceTerms/polarR6_ZoniGyro_CircularGeometry.h"

PolarR6_ZoniGyro_CircularGeometry::PolarR6_ZoniGyro_CircularGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid)
    , Rmax(Rmax)
{
}

double PolarR6_ZoniGyro_CircularGeometry::operator()(int i_r, int i_theta) const
{
    double r     = grid_.radius(i_r);
    double theta = grid_.theta(i_theta);
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * exp(tanh(10.0 * (r / Rmax) - 5.0)) *
               cos(11.0 * theta) -
           pow((r / Rmax), 4.0) *
               ((r / Rmax) *
                    (12.288 * (r / Rmax) * pow(((r / Rmax) - 1.0), 4.0) * cos(11.0 * theta) +
                     17.2032 * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta)) *
                    exp(-tanh(10.0 * (r / Rmax) - 5.0)) +
                (r / Rmax) *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                    (10.0 * pow(tanh(10.0 * (r / Rmax) - 5.0), 2.0) - 10.0) * exp(-tanh(10.0 * (r / Rmax) - 5.0)) -
                49.5616 * pow(((r / Rmax) - 1.0), 6.0) * exp(-tanh(10.0 * (r / Rmax) - 5.0)) * cos(11.0 * theta) +
                6.0 *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                    exp(-tanh(10.0 * (r / Rmax) - 5.0)));
}
