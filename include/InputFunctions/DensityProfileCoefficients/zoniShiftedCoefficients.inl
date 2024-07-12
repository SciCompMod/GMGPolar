#pragma once

#include "zoniShiftedCoefficients.h"

inline double ZoniShiftedCoefficients::alpha(const double& r) const {
    return exp(-tanh(20.0 * (r/Rmax) - 14.0));
}

inline double ZoniShiftedCoefficients::beta(const double& r) const {
    return 0.0;
}

inline double ZoniShiftedCoefficients::getAlphaJump() const {
    return alpha_jump;
}

