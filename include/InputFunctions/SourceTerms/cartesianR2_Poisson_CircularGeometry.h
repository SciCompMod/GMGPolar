#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR2_Poisson_CircularGeometry : public SourceTerm
{
public:
    CartesianR2_Poisson_CircularGeometry() = default;
    explicit CartesianR2_Poisson_CircularGeometry(double Rmax);
    virtual ~CartesianR2_Poisson_CircularGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
