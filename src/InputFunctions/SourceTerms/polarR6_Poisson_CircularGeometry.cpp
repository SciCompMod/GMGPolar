#include "../include/InputFunctions/SourceTerms/polarR6_Poisson_CircularGeometry.h"

PolarR6_Poisson_CircularGeometry::PolarR6_Poisson_CircularGeometry(const double& Rmax) : 
    Rmax(Rmax) 
{}

double PolarR6_Poisson_CircularGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const 
{
    return 0.0;
}