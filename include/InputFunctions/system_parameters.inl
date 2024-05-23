#pragma once

#include "system_parameters.h"

// In earlier versions denoted by 'rho_glob'
inline double rhs_f_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    static constexpr double Rmax = 1.3;
    static constexpr double map1_kappa = 0.3;
    static constexpr double map1_delta = 0.2;
    return (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) / (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) - ((r/Rmax) * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * ((-8.0) * M_PI * (r/Rmax) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 8.0 * M_PI * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta - 4.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * pow(sin_theta, 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 8.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) - 4.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * pow(cos_theta, 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta)) - 5.03290747193186 * (r/Rmax) * ((-2.0) * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta) / (208.641975308642 * pow(((r/Rmax) - 0.769230769230769), 2.0) + 1.0) + (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * ((-2.0) * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta) + (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * ((-4.0) * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * pow(sin_theta, 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 8.0 * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) - 4.0 * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * pow(cos_theta, 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta)) / (r/Rmax);
}

inline double u_D_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return 0.0;
}


// ------------- //
// Sonnendrucker //
// ------------- //

inline double alpha_Functor::operator()(const double& r) const {
    static constexpr double Rmax = 1.3;
    return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111);
}

inline double beta_Functor::operator()(const double& r) const {
    static constexpr double Rmax = 1.3;
    return pow((0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)), (double)((-1)));
    // return 0.0;
}

// ---- //
// Zoni //
// ---- //

// inline double alpha_Functor::operator()(const double& r) const {
//     static constexpr double Rmax = 1.3;
//     return exp(-tanh(10.0 * (r/Rmax) - 5.0));
// }

// inline double beta_Functor::operator()(const double& r) const {
//     static constexpr double Rmax = 1.3;
//     return exp(tanh(10.0 * (r/Rmax) - 5.0));
//     // return 0.0;
// }

// ------------ //
// Zoni-Shifted //
// ------------ //

// inline double alpha_Functor::operator()(const double& r) const {
//     static constexpr double Rmax = 1.3;
//     return exp(-tanh(20.0 * (r/Rmax) - 14.0));
// }

// inline double beta_Functor::operator()(const double& r) const {
//     static constexpr double Rmax = 1.3;
//     return exp(-tanh(20.0 * (r/Rmax) - 14.0));
//     // return 0.0;
// }

// ------- //
// Poisson // 
// ------- //

// inline double alpha_Functor::operator()(const double& r) const {
//     return 0.0;
// }

// inline double beta_Functor::operator()(const double& r) const {
//     return 0.0;
// }

