#pragma once

#include "sonnendruckerGyroCoefficients.h"

inline double SonnendruckerGyroCoefficients::alpha(const double& r) const {
    return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111);
}

inline double SonnendruckerGyroCoefficients::beta(const double& r) const {
    return pow((0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)), (double)((-1)));
}

inline double SonnendruckerGyroCoefficients::getAlphaJump() const {
    return alpha_jump;
}


