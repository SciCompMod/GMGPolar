#include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"

CzarnyGeometry::CzarnyGeometry()
{
    initializeGeometry();
}

CzarnyGeometry::CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon,
                               const double& ellipticity_e)
    : Rmax(Rmax)
    , inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon)
    , ellipticity_e(ellipticity_e)
{
    initializeGeometry();
}

void CzarnyGeometry::initializeGeometry()
{
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}
