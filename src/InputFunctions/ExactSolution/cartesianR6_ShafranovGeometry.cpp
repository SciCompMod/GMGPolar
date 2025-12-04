#include "../include/InputFunctions/ExactSolution/cartesianR6_ShafranovGeometry.h"

CartesianR6_ShafranovGeometry::CartesianR6_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta)
    : Rmax(Rmax)
    , elongation_kappa(elongation_kappa)
    , shift_delta(shift_delta)
{
}

double CartesianR6_ShafranovGeometry::exact_solution(double r, double theta) const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return 0.4096 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
           sin(M_PI * (2.0 * elongation_kappa * (r / Rmax) * sin_theta + 2.0 * (r / Rmax) * sin_theta)) *
           cos(M_PI * ((-2.0) * shift_delta * ((r / Rmax) * (r / Rmax)) -
                       2.0 * elongation_kappa * (r / Rmax) * cos_theta + 2.0 * (r / Rmax) * cos_theta));
}
