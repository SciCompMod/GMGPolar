#pragma once

#include <cmath>

#include "../exactSolution.h"

class Refined_CircularGeometry : public ExactSolution { 
public:
    Refined_CircularGeometry() = default;
    explicit Refined_CircularGeometry(const double& Rmax);
    virtual ~Refined_CircularGeometry() = default;

    double exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};
