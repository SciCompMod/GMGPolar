#pragma once

#include <cmath>

#include "../exactSolution.h"

class PolarR6_CulhamGeometry : public ExactSolution
{
public:
    PolarR6_CulhamGeometry() = default;
    explicit PolarR6_CulhamGeometry(double Rmax);
    virtual ~PolarR6_CulhamGeometry() = default;

    double exact_solution(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
