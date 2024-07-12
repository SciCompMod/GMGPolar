#pragma once

#include "zoniCoefficients.h"

inline double ZoniCoefficients::alpha(const double& r) const {
    return exp(-tanh(10.0 * (r/Rmax) - 5.0));
}

inline double ZoniCoefficients::beta(const double& r) const {
    return 0.0;
}

inline double ZoniCoefficients::getAlphaJump() const {
    return alpha_jump;
}

