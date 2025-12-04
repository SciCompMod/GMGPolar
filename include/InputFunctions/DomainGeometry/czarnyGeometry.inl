#pragma once

#include "czarnyGeometry.h"

// In earlier versions denoted by 'x'
inline double CzarnyGeometry::Fx(double r, double theta) const {
    double cos_theta = std::cos(theta);
    return (1.0 - sqrt(1.0 + inverse_aspect_ratio_epsilon * (inverse_aspect_ratio_epsilon + 2.0 * (r/Rmax) * cos_theta))) / inverse_aspect_ratio_epsilon;
}

// In earlier versions denoted by 'y'
inline double CzarnyGeometry::Fy(double r, double theta) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return ellipticity_e * factor_xi * (r/Rmax) * sin_theta / (1.0 + (1.0 - sqrt(1.0 + inverse_aspect_ratio_epsilon * (inverse_aspect_ratio_epsilon + 2.0 * (r/Rmax) * cos_theta))));
}

// In earlier versions denoted by 'Jrr'
inline double CzarnyGeometry::dFx_dr(double r, double theta) const {
    double cos_theta = std::cos(theta);
    return - (cos_theta) / (Rmax * sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r/Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0));
}

// In earlier versions denoted by 'Jtr'
inline double CzarnyGeometry::dFy_dr(double r, double theta) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double temp = sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r/Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0);
    return (ellipticity_e * factor_xi * sin_theta) / (Rmax * (2.0 - temp)) + (ellipticity_e * factor_xi * inverse_aspect_ratio_epsilon * r * sin_theta * cos_theta) / (Rmax * Rmax * temp * (2.0 - temp) * (2.0 - temp));
}

// In earlier versions denoted by 'Jrt'
inline double CzarnyGeometry::dFx_dt(double r, double theta) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return ((r/Rmax) * sin_theta) / (sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r/Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0));
}

// In earlier versions denoted by 'Jtt'
inline double CzarnyGeometry::dFy_dt(double r, double theta) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double temp = sqrt(inverse_aspect_ratio_epsilon * (2.0 * (r/Rmax) * cos_theta + inverse_aspect_ratio_epsilon) + 1.0);
    return (ellipticity_e * factor_xi * (r/Rmax) * cos_theta) / (2.0 - temp) - (ellipticity_e * factor_xi * inverse_aspect_ratio_epsilon * (r/Rmax) * (r/Rmax) * sin_theta * sin_theta) / (temp * (2.0 - temp) * (2.0 - temp));
}
