#pragma once

#include <cmath>

#include "../exactSolution.h"

class CartesianR6_CircularGeometry : public ExactSolution
{
public:
    CartesianR6_CircularGeometry() = default;
    explicit CartesianR6_CircularGeometry(double Rmax);
    virtual ~CartesianR6_CircularGeometry() = default;

    double exact_solution(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
