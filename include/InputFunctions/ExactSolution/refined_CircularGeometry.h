#pragma once

#include <cmath>

#include "../exactSolution.h"

class Refined_CircularGeometry : public ExactSolution
{
public:
    Refined_CircularGeometry() = default;
    explicit Refined_CircularGeometry(double Rmax);
    virtual ~Refined_CircularGeometry() = default;

    double exact_solution(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
