#pragma once

#include <cmath>

#include "../exactSolution.h"

class PolarR6_CircularGeometry : public ExactSolution
{
public:
    PolarR6_CircularGeometry() = default;
    explicit PolarR6_CircularGeometry(const double& Rmax);
    virtual ~PolarR6_CircularGeometry() = default;

    double exact_solution(const double& r, const double& theta, const double& sin_theta,
                          const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};
