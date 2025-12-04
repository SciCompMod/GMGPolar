#pragma once

#include <cmath>

#include "../exactSolution.h"

class Refined_CulhamGeometry : public ExactSolution
{
public:
    Refined_CulhamGeometry() = default;
    explicit Refined_CulhamGeometry(double Rmax);
    virtual ~Refined_CulhamGeometry() = default;

    double exact_solution(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
