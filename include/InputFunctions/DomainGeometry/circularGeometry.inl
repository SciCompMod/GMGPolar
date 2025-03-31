#pragma once

#include "circularGeometry.h"

// In earlier versions denoted by 'x'
inline double CircularGeometry::Fx(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (r/Rmax) * cos_theta;
}

// In earlier versions denoted by 'y'
inline double CircularGeometry::Fy(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (r/Rmax) * sin_theta;
}


// In earlier versions denoted by 'Jrr'
inline double CircularGeometry::dFx_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (cos_theta) / Rmax;
}

// In earlier versions denoted by 'Jtr'
inline double CircularGeometry::dFy_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (sin_theta) / Rmax;
}

// In earlier versions denoted by 'Jrt'
inline double CircularGeometry::dFx_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (-(r/Rmax)) * sin_theta;
}

// In earlier versions denoted by 'Jtt'
inline double CircularGeometry::dFy_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (r/Rmax) * cos_theta;
}
