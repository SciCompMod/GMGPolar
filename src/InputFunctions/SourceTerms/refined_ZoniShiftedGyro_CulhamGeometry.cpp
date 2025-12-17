#include "../include/InputFunctions/SourceTerms/refined_ZoniShiftedGyro_CulhamGeometry.h"

Refined_ZoniShiftedGyro_CulhamGeometry::Refined_ZoniShiftedGyro_CulhamGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid)
    , Rmax(Rmax)
{
}

double Refined_ZoniShiftedGyro_CulhamGeometry::operator()(int i_r, int i_theta) const
{
    double r         = grid_.radius(i_r);
    double theta     = grid_.theta(i_theta);
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) - 0.0 * (r / Rmax) - 0.0 +
            exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
               cos(21.0 * theta) +
           (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) - 4.00652973929511e-05 +
            exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
               cos(9.0 * theta);
}