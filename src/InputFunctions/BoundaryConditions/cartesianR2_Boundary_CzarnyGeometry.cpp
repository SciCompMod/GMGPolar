#include "../include/InputFunctions/BoundaryConditions/cartesianR2_Boundary_CzarnyGeometry.h"

void CartesianR2_Boundary_CzarnyGeometry::initializeGeometry(){
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}

CartesianR2_Boundary_CzarnyGeometry::CartesianR2_Boundary_CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e) : 
    Rmax(Rmax),
    inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon),
    ellipticity_e(ellipticity_e)
{
    initializeGeometry();
}

double CartesianR2_Boundary_CzarnyGeometry::u_D(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    double temp = sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r/Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0);
    return (1.0 - (r/Rmax) * (r/Rmax)) * 
        sin(2.0 * M_PI * ellipticity_e * (r/Rmax) * sin_theta * factor_xi / (2.0 - temp)) * 
        cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon);
}


double CartesianR2_Boundary_CzarnyGeometry::u_D_Interior(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    double temp = sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r/Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0);
    return (1.0 - (r/Rmax) * (r/Rmax)) * 
        sin(2.0 * M_PI * ellipticity_e * (r/Rmax) * sin_theta * factor_xi / (2.0 - temp)) * 
        cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon);
}
