#pragma once

#include <cmath>

#include "../exactSolution.h"

namespace gmgpolar
{

class CartesianR6_ShafranovGeometry
{
public:
    CartesianR6_ShafranovGeometry() = default;
    explicit CartesianR6_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);
    virtual ~CartesianR6_ShafranovGeometry() = default;

    double exact_solution(double r, double theta) const;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
} // namespace gmgpolar
