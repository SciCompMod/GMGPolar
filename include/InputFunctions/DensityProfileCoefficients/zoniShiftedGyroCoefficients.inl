#pragma once

#include "zoniShiftedGyroCoefficients.h"

inline double ZoniShiftedGyroCoefficients::alpha(const double& r) const {
    return exp(-tanh(20.0 * (r/Rmax) - 14.0));
}

inline double ZoniShiftedGyroCoefficients::beta(const double& r) const {
    return exp(tanh(20.0 * (r/Rmax) - 14.0));
}

inline double ZoniShiftedGyroCoefficients::getAlphaJump() const {
    return alpha_jump;
}


