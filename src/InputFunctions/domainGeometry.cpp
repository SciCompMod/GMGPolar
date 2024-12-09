#include "../../include/InputFunctions/domainGeometry.h"

DomainGeometry::DomainGeometry()
{
    initializeGeometry();
}

DomainGeometry::DomainGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e)
    : Rmax(Rmax), inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon), ellipticity_e(ellipticity_e) {
    initializeGeometry();
}

void DomainGeometry::initializeGeometry() {
    assert(inverse_aspect_ratio_epsilon < 2.0); // Ensure the epsilon is in a valid range
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}

