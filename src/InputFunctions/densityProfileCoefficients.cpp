#include "../../include/InputFunctions/densityProfileCoefficients.h"

DensityProfileCoefficients::DensityProfileCoefficients(const double& Rmax, const double& alpha_jump) : 
    Rmax(Rmax),
    alpha_jump(alpha_jump) 
{}

double DensityProfileCoefficients::getAlphaJump() const {
    return alpha_jump;
}