#include "../include/InputFunctions/ExactSolution/refined_CulhamGeometry.h"

Refined_CulhamGeometry::Refined_CulhamGeometry(const double& Rmax) : 
    Rmax(Rmax) 
{}

double Refined_CulhamGeometry::exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return 0.0;
}
