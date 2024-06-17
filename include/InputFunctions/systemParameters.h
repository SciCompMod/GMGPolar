#pragma once

#include <cmath>

class SystemParameters {
public:
    explicit SystemParameters() = default;

    double alpha(const double& r) const;
    double beta(const double& r) const;
    double getAlphaJump() const;

    double u_D(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const;
    double u_D_Interior(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const;

    double rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const;

private:
    const double Rmax = 1.3;

    /* Shafranov Geometry */
    // const double map1_kappa = 0.3;
    // const double map1_delta = 0.2;

    /* Czarny/Triangular Geometry */
    const double map2_epsilon = 0.3;
    const double map2_e = 1.4;

    // const double alpha_jump = 0.5; /* Poisson */
    const double alpha_jump = 0.66; /* Sonnendrucker */
    // const double alpha_jump = 0.4837; /* Zoni */
    // const double alpha_jump = 0.7081; /* Zoni Shifted */
};

#include "systemParameters.inl"