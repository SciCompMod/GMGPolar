#pragma once

#include <cmath>

#include "../exactSolution.h"

class Refined_CulhamGeometry : public ExactSolution
{
public:
    Refined_CulhamGeometry() = default;
    explicit Refined_CulhamGeometry(const double& Rmax);
    virtual ~Refined_CulhamGeometry() = default;

    double exact_solution(const double& r, const double& theta, const double& sin_theta,
                          const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};
