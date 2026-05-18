#pragma once

#include <cmath>

#include "../exactSolution.h"

namespace gmgpolar
{

class CartesianR6_CzarnyGeometry
{
public:
    explicit CartesianR6_CzarnyGeometry();
    explicit CartesianR6_CzarnyGeometry(double Rmax, double inverse_aspect_ratio_epsilon, double ellipticity_e);

    virtual ~CartesianR6_CzarnyGeometry() = default;

    double exact_solution(double r, double theta) const;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
} // namespace gmgpolar
