#include "../../include/Operator/operator.h"

Operator::Operator(const GMGPolar& gmgpolar, const PolarGrid& grid) : 
    kappa_eps_(gmgpolar.kappa_eps),
    delta_e_(gmgpolar.delta_e),
    Rmax_(gmgpolar.Rmax),
    geometry_(gmgpolar.geometry),
    sin_theta_(grid.ntheta()), 
    cos_theta_(grid.ntheta()) 
{
    
    #pragma omp parallel for
    for (int i = 0; i < grid.ntheta(); i++) {
        sin_theta_[i] = sin(grid.theta(i));
        cos_theta_[i] = cos(grid.theta(i));
    }
}