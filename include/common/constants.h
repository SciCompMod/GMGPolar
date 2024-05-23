#pragma once

enum alpha_coeff
{
    SONNENDRUCKER = 0,
    ZONI          = 1,
    ZONI_SHIFTED  = 2,
    POISSON       = 3
};

enum beta_coeff
{
    ZERO = 0,
    ALPHA_INVERSE = 1
};

enum geometry_type
{
    CIRCULAR   = 0,
    SHAFRANOV  = 1,
    TRIANGULAR = 2,
    CULHAM     = 3 
};

enum problem_type
{
    CARTESIAN_R2   = 0, 
    POLAR_R6       = 1,
    CARTESIAN_R6   = 2,
    REFINED_RADIUS = 3
};



// constexpr inline double Rmax = 1.3; 

// inline double F_x(const double& r, const double& theta, const double& sin_theta, const double& cos_theta)
// {
//     return (r/Rmax) * cos_theta;
// }

// inline double F_y(const double& r, const double& theta, const double& sin_theta, const double& cos_theta)
// {
//     return (r/Rmax) * sin_theta;
// }

//  double dFx_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta)
// {
//     return (cos_theta)/Rmax;
// }

//  double dFy_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta)
// {
//     return (sin_theta)/Rmax;
// }

//  double dFx_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta)
// {
//     return (-(r/Rmax)) * sin_theta;
// }

//  double dFy_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta)
// {
//     return (r/Rmax) * cos_theta;
// }

