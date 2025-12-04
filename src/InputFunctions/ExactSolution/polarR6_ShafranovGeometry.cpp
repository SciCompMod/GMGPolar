#include "../include/InputFunctions/ExactSolution/polarR6_ShafranovGeometry.h"

PolarR6_ShafranovGeometry::PolarR6_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta)
    : Rmax(Rmax)
    , elongation_kappa(elongation_kappa)
    , shift_delta(shift_delta)
{
}

double PolarR6_ShafranovGeometry::exact_solution(double r, double theta) const
{
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
