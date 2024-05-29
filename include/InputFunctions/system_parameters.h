#pragma once

#include <cmath>

struct alpha_Functor {
    double operator()(const double& r) const;
};

struct beta_Functor {
    double operator()(const double& r) const;
};

struct rhs_f_Functor {
    double operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const;
};

struct u_D_Functor {
    double operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const;
};

struct u_D_Interior_Functor {
    double operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const;
};

#include "system_parameters.inl"