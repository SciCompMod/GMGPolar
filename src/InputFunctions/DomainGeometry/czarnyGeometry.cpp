#include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"

CzarnyGeometry::CzarnyGeometry(){
    initializeGeometry();
}

CzarnyGeometry::CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e) : 
    Rmax(Rmax),
    inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon),
    ellipticity_e(ellipticity_e)
{
    initializeGeometry();
}

void CzarnyGeometry::initializeGeometry() {
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}

// // In earlier versions denoted by 'x'
// double CzarnyGeometry::Fx(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     return (1.0 - sqrt(1.0 + inverse_aspect_ratio_epsilon * (inverse_aspect_ratio_epsilon + 2.0 * (r/Rmax) * cos_theta))) / inverse_aspect_ratio_epsilon;
// }

// // In earlier versions denoted by 'y'
// double CzarnyGeometry::Fy(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     return ellipticity_e * factor_xi * (r/Rmax) * sin_theta / (1.0 + (1.0 - sqrt(1.0 + inverse_aspect_ratio_epsilon * (inverse_aspect_ratio_epsilon + 2.0 * (r/Rmax) * cos_theta))));
// }

// // In earlier versions denoted by 'Jrr'
// double CzarnyGeometry::dFx_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     return - (cos_theta) / (Rmax * sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r/Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0));
// }

// // In earlier versions denoted by 'Jtr'
// double CzarnyGeometry::dFy_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     double temp = sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r/Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0);
//     return (ellipticity_e * factor_xi * sin_theta) / (Rmax * (2.0 - temp)) + (ellipticity_e * factor_xi * inverse_aspect_ratio_epsilon * r * sin_theta * cos_theta) / (Rmax * Rmax * temp * (2.0 - temp) * (2.0 - temp));
// }

// // In earlier versions denoted by 'Jrt'
// double CzarnyGeometry::dFx_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     return ((r/Rmax) * sin_theta) / (sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r/Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0));
// }

// // In earlier versions denoted by 'Jtt'
// double CzarnyGeometry::dFy_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     double temp = sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r/Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0);
//     return (ellipticity_e * factor_xi * (r/Rmax) * cos_theta) / (2.0 - temp) - (ellipticity_e * factor_xi * inverse_aspect_ratio_epsilon * (r/Rmax) * (r/Rmax) * sin_theta * sin_theta) / (temp * (2.0 - temp) * (2.0 - temp));
// }
