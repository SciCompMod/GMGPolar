#include "../include/InputFunctions/DensityProfileCoefficients/sonnendruckerGyroCoefficients.h"

#include <iostream>

SonnendruckerGyroCoefficients::SonnendruckerGyroCoefficients(const double& _Rmax, const double& _alpha_jump) : 
    Rmax(_Rmax),
    alpha_jump(_alpha_jump) 
{}
