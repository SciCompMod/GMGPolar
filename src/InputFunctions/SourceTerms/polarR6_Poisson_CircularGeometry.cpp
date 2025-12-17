#include "../include/InputFunctions/SourceTerms/polarR6_Poisson_CircularGeometry.h"

PolarR6_Poisson_CircularGeometry::PolarR6_Poisson_CircularGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid) , Rmax(Rmax)
{
}

double PolarR6_Poisson_CircularGeometry::operator()(int i_r, int i_theta) const
{
    return (-pow((r / Rmax), 4.0)) * (14.7456 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                                      1.0 * (r / Rmax) *
                                          (12.288 * (r / Rmax) * pow(((r / Rmax) - 1.0), 4.0) * cos(11.0 * theta) +
                                           17.2032 * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta)) -
                                      34.816 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta));
}
