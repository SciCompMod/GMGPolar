#include "../../include/InputFunctions/boundaryConditions.h"

#ifdef GEOM_SHAFRANOV   
BoundaryConditions::BoundaryConditions()
{
}

BoundaryConditions::BoundaryConditions(const double& Rmax, const double& elongation_kappa, const double& shift_delta)
    : Rmax(Rmax), elongation_kappa(elongation_kappa), shift_delta(shift_delta)
{
}    
#else
BoundaryConditions::BoundaryConditions()
{
    initializeGeometry();
}


BoundaryConditions::BoundaryConditions(const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e)
    : Rmax(Rmax), inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon), ellipticity_e(ellipticity_e) {
    initializeGeometry();
}

void BoundaryConditions::initializeGeometry() {
    assert(inverse_aspect_ratio_epsilon < 2.0); // Ensure epsilon is in a valid range
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}
#endif
