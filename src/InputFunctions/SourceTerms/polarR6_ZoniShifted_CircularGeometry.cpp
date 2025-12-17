#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShifted_CircularGeometry.h"

PolarR6_ZoniShifted_CircularGeometry::PolarR6_ZoniShifted_CircularGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid) , Rmax(Rmax)
{
}

double PolarR6_ZoniShifted_CircularGeometry::operator()(int i_r, int i_theta) const
{
    return (-pow((r / Rmax), 4.0)) *
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
