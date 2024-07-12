#pragma once

#include "poissonCoefficients.h"

inline double PoissonCoefficients::alpha(const double& r) const {
    return 1.0;
}

inline double PoissonCoefficients::beta(const double& r) const {
    return 0.0;
}

inline double PoissonCoefficients::getAlphaJump() const {
    return alpha_jump;
}
