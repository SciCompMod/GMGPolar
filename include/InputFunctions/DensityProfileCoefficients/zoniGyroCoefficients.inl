#pragma once

#include "zoniGyroCoefficients.h"

inline double ZoniGyroCoefficients::alpha(const double& r) const {
    return exp(-tanh(10.0 * (r/Rmax) - 5.0));
}

inline double ZoniGyroCoefficients::beta(const double& r) const {
    return exp(tanh(10.0 * (r/Rmax) - 5.0));
}

inline double ZoniGyroCoefficients::getAlphaJump() const {
    return alpha_jump;
}


