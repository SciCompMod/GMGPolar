#include "../include/InputFunctions/DensityProfileCoefficients/sonnendruckerGyroCoefficients.h"

SonnendruckerGyroCoefficients::SonnendruckerGyroCoefficients(const double& _Rmax, const double& _alpha_jump) : 
    Rmax(_Rmax),
    alpha_jump(_alpha_jump) 
{}

double SonnendruckerGyroCoefficients::alpha(const double& r) const {
    return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111);
}

double SonnendruckerGyroCoefficients::beta(const double& r) const {
    return pow((0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)), (double)((-1)));
}

double SonnendruckerGyroCoefficients::getAlphaJump() const {
    return alpha_jump;
}

