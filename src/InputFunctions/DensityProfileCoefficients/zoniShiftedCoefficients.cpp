#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h"
#include <iostream>
ZoniShiftedCoefficients::ZoniShiftedCoefficients(const double& Rmax, const double& alpha_jump) : 
    Rmax(Rmax),
    alpha_jump(alpha_jump) 
{}

double ZoniShiftedCoefficients::alpha(const double& r) const {
    return exp(-tanh(20.0 * (r/Rmax) - 14.0));
}

double ZoniShiftedCoefficients::beta(const double& r) const {
    return 0.0;
}

double ZoniShiftedCoefficients::getAlphaJump() const {
    return alpha_jump;
}

