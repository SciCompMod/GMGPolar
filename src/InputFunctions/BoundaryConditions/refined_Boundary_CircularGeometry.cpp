#include "../include/InputFunctions/BoundaryConditions/refined_Boundary_CircularGeometry.h"

Refined_Boundary_CircularGeometry::Refined_Boundary_CircularGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double Refined_Boundary_CircularGeometry::u_D(const double& r, const double& theta)const
{
    return ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) - 0.0 * (r / Rmax) - 0.0 +
            exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
               cos(21.0 * theta) +
           (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) - 4.00652973929511e-05 +
            exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
               cos(9.0 * theta);
}

double Refined_Boundary_CircularGeometry::u_D_Interior(const double& r, const double& theta)const
{
    return ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) - 0.0 * (r / Rmax) - 0.0 +
            exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
               cos(21.0 * theta) +
           (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) - 4.00652973929511e-05 +
            exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
               cos(9.0 * theta);
}
