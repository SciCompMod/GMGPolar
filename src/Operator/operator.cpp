#include "../../include/Operator/operator.h"

Operator::Operator(const GMGPolar& gmgpolar, const PolarGrid& grid) : 
    kappa_eps_(gmgpolar.kappa_eps),
    delta_e_(gmgpolar.delta_e),
    Rmax_(gmgpolar.Rmax),
    geometry_(gmgpolar.geometry),
    sin_theta_(grid.ntheta()), 
    cos_theta_(grid.ntheta()),
    dFx_dr_(gmgpolar.dFx_dr_),
    dFy_dr_(gmgpolar.dFy_dr_),
    dFx_dt_(gmgpolar.dFx_dt_),
    dFy_dt_(gmgpolar.dFy_dt_),
    alpha_(gmgpolar.alpha_),
    beta_(gmgpolar.beta_),

    numberMod0Circles((grid.numberSmootherCircles() + 2) / 3),
    numberMod1Circles((grid.numberSmootherCircles() + 1) / 3),
    numberMod2Circles((grid.numberSmootherCircles() + 0) / 3),

    numberDiv3Radials((grid.ntheta() + 0) / 3),

    Mod1CircleMutexes(numberMod1Circles),
    Mod1CircleDepCounter(numberMod1Circles),

    Mod2CircleMutexes(numberMod2Circles),
    Mod2CircleDepCounter(numberMod2Circles),

    Mod1RadialMutexes(numberDiv3Radials),
    Mod1RadialDepCounter(numberDiv3Radials),

    Mod2RadialMutexes(numberDiv3Radials),
    Mod2RadialDepCounter(numberDiv3Radials),

    additional_circle_task(grid.numberSmootherCircles() % 3),
    additional_radial_task(grid.ntheta() % 3)

{
    #pragma omp parallel for
    for (int i = 0; i < grid.ntheta(); i++) {
        sin_theta_[i] = sin(grid.theta(i));
        cos_theta_[i] = cos(grid.theta(i));
    }
}