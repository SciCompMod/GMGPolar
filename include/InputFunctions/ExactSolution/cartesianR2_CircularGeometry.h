#pragma once

#include <cmath>

#include "../exactSolution.h"

namespace gmgpolar
{

class CartesianR2_CircularGeometry : public ExactSolution
{
public:
    CartesianR2_CircularGeometry() = default;
    explicit CartesianR2_CircularGeometry(double Rmax);
    virtual ~CartesianR2_CircularGeometry() = default;

    double exact_solution(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
} // namespace gmgpolar
