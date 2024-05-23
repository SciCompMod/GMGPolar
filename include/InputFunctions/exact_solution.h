#pragma once

#include <cmath>

struct exact_solution_Functor { // In earlier versions denoted by 'exact_phi'
    double operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const;
};

#include "exact_solution.inl"