#pragma once

#include <cmath>

class ExactSolution {
public:
    explicit ExactSolution() = default;

    double exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const;
private:
    const double Rmax = 1.3;

    /* Shafranov Geometry */
    const double map1_kappa = 0.3;
    const double map1_delta = 0.2;

    /* Czarny/Triangular Geometry */
    // const double map2_epsilon = 0.3;
    // const double map2_e = 1.4;
};

#include "exactSolution.inl"