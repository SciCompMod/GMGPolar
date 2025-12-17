#include "../include/InputFunctions/SourceTerms/cartesianR2_Poisson_CircularGeometry.h"

CartesianR2_Poisson_CircularGeometry::CartesianR2_Poisson_CircularGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid)
    , Rmax(Rmax)
{
}

double CartesianR2_Poisson_CircularGeometry::operator()(int i_r, int i_theta) const
{
    double r         = grid_.radius(i_r);
    double theta     = grid_.theta(i_theta);
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return 8.0 * M_PI * (r / Rmax) * sin_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) *
               cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
           8.0 * M_PI * (r / Rmax) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
               sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
           8.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * pow(sin_theta, 2.0) *
               sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
           8.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
               pow(cos_theta, 2.0) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
           4.0 * sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta);
}
