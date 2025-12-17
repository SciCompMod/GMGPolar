#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CulhamGeometry.h"

PolarR6_ZoniShiftedGyro_CulhamGeometry::PolarR6_ZoniShiftedGyro_CulhamGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid)
    , Rmax(Rmax)
{
}

double PolarR6_ZoniShiftedGyro_CulhamGeometry::operator()(int i_r, int i_theta) const
{
    double r         = grid_.radius(i_r);
    double theta     = grid_.theta(i_theta);
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}