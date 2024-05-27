#pragma once

#include "domain_geometry.h"

// ------------------ //
// Geometry: Circular //
// ------------------ //

// // In earlier versions denoted by 'x'
// inline double Fx_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3; 
//     return (r / Rmax) * cos_theta;
// }

// // In earlier versions denoted by 'y'
// inline double Fy_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3; 
//     return (r/Rmax) * sin_theta;
// }


// // In earlier versions denoted by 'Jrr'
// inline double dFx_dr_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     return (cos_theta)/Rmax;
// }

// // In earlier versions denoted by 'Jtr'
// inline double dFy_dr_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     return (sin_theta)/Rmax;
// }

// // In earlier versions denoted by 'Jrt'
// inline double dFx_dt_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     return (-(r/Rmax)) * sin_theta;
// }

// // In earlier versions denoted by 'Jtt'
// inline double dFy_dt_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     return (r/Rmax) * cos_theta;
// }


// ------------------- //
// Geometry: Shafranov //
// ------------------- //

// // In earlier versions denoted by 'x'
// inline double Fx_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3; 
//     static constexpr double map1_kappa = 0.3;
//     static constexpr double map1_delta = 0.2;
//     return (-map1_delta) * ((r/Rmax) * (r/Rmax)) - map1_kappa * (r/Rmax) * cos_theta + (r/Rmax) * cos_theta;
// }

// // In earlier versions denoted by 'y'
// inline double Fy_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3; 
//     static constexpr double map1_kappa = 0.3;
//     static constexpr double map1_delta = 0.2;
//     return map1_kappa * (r/Rmax) * sin_theta + (r/Rmax) * sin_theta;
// }


// // In earlier versions denoted by 'Jrr'
// inline double dFx_dr_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     static constexpr double map1_kappa = 0.3;
//     static constexpr double map1_delta = 0.2;
//     return ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta)/Rmax;
// }

// // In earlier versions denoted by 'Jtr'
// inline double dFy_dr_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     static constexpr double map1_kappa = 0.3;
//     static constexpr double map1_delta = 0.2;
//     return ((map1_kappa + 1.0) * sin_theta)/Rmax;
// }

// // In earlier versions denoted by 'Jrt'
// inline double dFx_dt_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     static constexpr double map1_kappa = 0.3;
//     static constexpr double map1_delta = 0.2;
//     return (r/Rmax) * (map1_kappa * sin_theta - sin_theta);
// }

// // In earlier versions denoted by 'Jtt'
// inline double dFy_dt_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     static constexpr double map1_kappa = 0.3;
//     static constexpr double map1_delta = 0.2;
//     return (r/Rmax) * (map1_kappa * cos_theta + cos_theta);
// }


// --------------------------- //
// Geometry: Czarny/Triangular //
// --------------------------- //

// In earlier versions denoted by 'x'
inline double Fx_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    static constexpr double Rmax = 1.3; 
    static constexpr double map2_epsilon = 0.3;
    static constexpr double map2_e = 1.4;
    return (1.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)) / map2_epsilon;
}

// In earlier versions denoted by 'y'
inline double Fy_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    static constexpr double Rmax = 1.3; 
    static constexpr double map2_epsilon = 0.3;
    static constexpr double map2_e = 1.4;
    return map2_e * (r/Rmax) * sin_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)));
}


// In earlier versions denoted by 'Jrr'
inline double dFx_dr_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    static constexpr double Rmax = 1.3;
    static constexpr double map2_epsilon = 0.3;
    static constexpr double map2_e = 1.4;
    return ((-cos_theta) / sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0))/Rmax;
}

// In earlier versions denoted by 'Jtr'
inline double dFy_dr_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    static constexpr double Rmax = 1.3;
    static constexpr double map2_epsilon = 0.3;
    static constexpr double map2_e = 1.4;
    return (map2_e * map2_epsilon * (r/Rmax) * sin_theta * cos_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * pow((2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)), 2.0) * sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)) + map2_e * sin_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0))))/Rmax;
}

// In earlier versions denoted by 'Jrt'
inline double dFx_dt_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    static constexpr double Rmax = 1.3;
    static constexpr double map2_epsilon = 0.3;
    static constexpr double map2_e = 1.4;
    return (r/Rmax) * sin_theta / sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0);
}

// In earlier versions denoted by 'Jtt'
inline double dFy_dt_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    static constexpr double Rmax = 1.3;
    static constexpr double map2_epsilon = 0.3;
    static constexpr double map2_e = 1.4;
    return (r/Rmax) * ((-map2_e) * map2_epsilon * (r/Rmax) * pow(sin_theta, 2.0) / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * pow((2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)), 2.0) * sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)) + map2_e * cos_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0))));
}

// ------ //
// Culham //
// ------ //

// inline double TransformationHelper::my_sum(std::array<double, 1001>& f, int64_t start_idx, int64_t end_idx) const
// {
//     int64_t i;
//     double result;
//     result = 0.0;
//     #pragma omp parallel for reduction(+: result)
//     for (i = start_idx; i < end_idx; i += 1)
//     {
//         result += f[i];
//     }
//     return result;
// }

// inline double TransformationHelper::q(double rr) const
// {
//     return 0.8 - 0.1 * (rr * rr);
// }

// inline double TransformationHelper::dq(double rr) const
// {
//     return (-0.2) * rr;
// }

// inline double TransformationHelper::p(double rr) const
// {
//     return 100000.0 - 90000.0 * (rr * rr);
// }

// inline double TransformationHelper::dp(double rr) const
// {
//     return (-180000.0) * rr;
// }

// inline double TransformationHelper::dg(double rr, double g) const
// {
//     return ((-g) * (0.0625000000000001 * (rr * rr) / pow((1.0 - 0.125 * (rr * rr)), 2.0) + 2.0 / (4.0 - 0.5 * (rr * rr))) + 2.261946711816e-06 * (4.0 - 0.5 * (rr * rr)) / g) / (rr / (4.0 - 0.5 * (rr * rr)) + (4.0 - 0.5 * (rr * rr)) / (g * rr));
// }

// inline double TransformationHelper::double_deriv(double rr, double c, double g, double dg, double val, double d_val) const
// {
//     return c * val / (rr * rr) - d_val * (pow(rr, (double)((-1))) + (4.0 - 0.5 * (rr * rr)) * (2.0 * dg * rr / (4.0 - 0.5 * (rr * rr)) + 0.125 * g * (rr * rr) / pow((1.0 - 0.125 * (rr * rr)), 2.0) + 2.0 * g / (4.0 - 0.5 * (rr * rr))) / (g * rr));
// }

// inline double TransformationHelper::g(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return g_array[ri];
//     }
//     else
//     {
//         m = (g_array[ri + 1] - g_array[ri]) / dr;
//         c = g_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// inline double TransformationHelper::Delta(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return Delta_array[ri];
//     }
//     else
//     {
//         m = (Delta_array[ri + 1] - Delta_array[ri]) / dr;
//         c = Delta_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// inline double TransformationHelper::Delta_prime(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return Delta_prime_array[ri];
//     }
//     else
//     {
//         m = (Delta_prime_array[ri + 1] - Delta_prime_array[ri]) / dr;
//         c = Delta_prime_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// inline double TransformationHelper::E(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return E_array[ri];
//     }
//     else
//     {
//         m = (E_array[ri + 1] - E_array[ri]) / dr;
//         c = E_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// inline double TransformationHelper::T(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return T_array[ri];
//     }
//     else
//     {
//         m = (T_array[ri + 1] - T_array[ri]) / dr;
//         c = T_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// inline double TransformationHelper::E_prime(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return E_prime_array[ri];
//     }
//     else
//     {
//         m = (E_prime_array[ri + 1] - E_prime_array[ri]) / dr;
//         c = E_prime_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// inline double TransformationHelper::T_prime(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return T_prime_array[ri];
//     }
//     else
//     {
//         m = (T_prime_array[ri + 1] - T_prime_array[ri]) / dr;
//         c = T_prime_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// inline double TransformationHelper::P(double rr) const
// {
//     if (rr == 0)
//     {
//         return 0.0;
//     }
//     else
//     {
//         return 0.005 * pow(rr, 3.0) + 0.1 * rr * Delta(rr) - 0.1 * pow(E(rr), 2.0) - pow(T(rr), 2.0) / rr;
//     }
// }

// inline double TransformationHelper::dP(double rr) const
// {
//     if (rr == 0)
//     {
//         return 0.0;
//     }
//     else
//     {
//         return 0.015 * (rr * rr) + 0.1 * rr * Delta_prime(rr) + 0.1 * Delta(rr) - 0.2 * E(rr) * E_prime(rr) - 2.0 * T(rr) * T_prime(rr) / rr + pow(T(rr), 2.0) / (rr * rr);
//     }
// }

// // In earlier versions denoted by 'x'
// inline double Fx_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3; 
//     return (r/Rmax) * cos_theta + helper_->Delta((r/Rmax)) - helper_->E((r/Rmax)) * cos_theta - helper_->P((r/Rmax)) * cos_theta + helper_->T((r/Rmax)) * cos(2.0 * theta) + 5.0;
// }

// // In earlier versions denoted by 'y'
// inline double Fy_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3; 
//     return (r/Rmax) * sin_theta - helper_->E((r/Rmax)) * sin_theta - helper_->P((r/Rmax)) * sin_theta - helper_->T((r/Rmax)) * sin(2.0 * theta);
// }


// // In earlier versions denoted by 'Jrr'
// inline double dFx_dr_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     return (helper_->Delta_prime((r/Rmax)) - helper_->E_prime((r/Rmax)) * cos_theta + helper_->T_prime((r/Rmax)) * cos(2.0 * theta) - helper_->dP((r/Rmax)) * cos_theta + cos_theta)/Rmax;
// }

// // In earlier versions denoted by 'Jtr'
// inline double dFy_dr_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     return ((-helper_->E_prime((r/Rmax))) * sin_theta - helper_->T_prime((r/Rmax)) * sin(2.0 * theta) - helper_->dP((r/Rmax)) * sin_theta + sin_theta)/Rmax;
// }

// // In earlier versions denoted by 'Jrt'
// inline double dFx_dt_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     return (-(r/Rmax)) * sin_theta + helper_->E((r/Rmax)) * sin_theta + helper_->P((r/Rmax)) * sin_theta - 2.0 * helper_->T((r/Rmax)) * sin(2.0 * theta);
// }

// // In earlier versions denoted by 'Jtt'
// inline double dFy_dt_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     static constexpr double Rmax = 1.3;
//     return (r/Rmax) * cos_theta - helper_->E((r/Rmax)) * cos_theta - helper_->P((r/Rmax)) * cos_theta - 2.0 * helper_->T((r/Rmax)) * cos(2.0 * theta);
// }