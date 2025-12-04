#pragma once

#include <cmath>

#include "../exactSolution.h"

class CartesianR2_CircularGeometry : public ExactSolution
{
public:
    CartesianR2_CircularGeometry() = default;
    explicit CartesianR2_CircularGeometry(const double& Rmax);
    virtual ~CartesianR2_CircularGeometry() = default;

    double exact_solution(const double& r, const double& theta) const override;

private:
    const double Rmax = 1.3;
};
