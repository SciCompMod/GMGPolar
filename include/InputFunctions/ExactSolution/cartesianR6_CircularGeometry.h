#pragma once

#include <cmath>

#include "../exactSolution.h"

class CartesianR6_CircularGeometry : public ExactSolution { 
public:
    CartesianR6_CircularGeometry() = default;
    explicit CartesianR6_CircularGeometry(const double& Rmax);
    virtual ~CartesianR6_CircularGeometry() = default;

    double exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};
