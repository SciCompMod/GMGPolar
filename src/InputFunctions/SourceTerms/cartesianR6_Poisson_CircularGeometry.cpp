#include "../include/InputFunctions/SourceTerms/cartesianR6_Poisson_CircularGeometry.h"

CartesianR6_Poisson_CircularGeometry::CartesianR6_Poisson_CircularGeometry(const double& Rmax) : 
    Rmax(Rmax) 
{}

double CartesianR6_Poisson_CircularGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const 
{
    return 0.0;
}