#pragma once

#include "circularGeometry.h"

// In earlier versions denoted by 'x'
inline double CircularGeometry::Fx(double r, double theta) const {
    return (r/Rmax) * std::cos(theta);
}

// In earlier versions denoted by 'y'
inline double CircularGeometry::Fy(double r, double theta) const {
    return (r/Rmax) * std::sin(theta);
}


// In earlier versions denoted by 'Jrr'
inline double CircularGeometry::dFx_dr(double r, double theta) const {
    return (std::cos(theta)) / Rmax;
}

// In earlier versions denoted by 'Jtr'
inline double CircularGeometry::dFy_dr(double r, double theta) const {
    return (std::sin(theta)) / Rmax;
}

// In earlier versions denoted by 'Jrt'
inline double CircularGeometry::dFx_dt(double r, double theta) const {
    return (-(r/Rmax)) * std::sin(theta);
}

// In earlier versions denoted by 'Jtt'
inline double CircularGeometry::dFy_dt(double r, double theta) const {
    return (r/Rmax) * std::cos(theta);
}
