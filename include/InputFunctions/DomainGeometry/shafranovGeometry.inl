#pragma once

#include "shafranovGeometry.h"


// In earlier versions denoted by 'x'
inline double ShafranovGeometry::Fx(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (1.0 - elongation_kappa) * (r/Rmax) * cos_theta - shift_delta * (r/Rmax) * (r/Rmax);
}

// In earlier versions denoted by 'y'
inline double ShafranovGeometry::Fy(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (1.0 + elongation_kappa) * (r/Rmax) * sin_theta;
}


// In earlier versions denoted by 'Jrr'
inline double ShafranovGeometry::dFx_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return ( (Rmax - elongation_kappa * Rmax) * cos_theta - 2.0 * shift_delta * r ) / ( Rmax * Rmax );
}

// In earlier versions denoted by 'Jtr'
inline double ShafranovGeometry::dFy_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (elongation_kappa + 1.0) * sin_theta / Rmax;
}

// In earlier versions denoted by 'Jrt'
inline double ShafranovGeometry::dFx_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return ((elongation_kappa - 1.0) * r * sin_theta) / Rmax;
}

// In earlier versions denoted by 'Jtt'
inline double ShafranovGeometry::dFy_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return ((elongation_kappa + 1.0) * r * cos_theta) / Rmax;
}