#pragma once

#include <cmath>

#include "../exactSolution.h"

class PolarR6_CircularGeometry : public ExactSolution
{
public:
    PolarR6_CircularGeometry() = default;
    explicit PolarR6_CircularGeometry(double Rmax);
    virtual ~PolarR6_CircularGeometry() = default;

    double exact_solution(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
