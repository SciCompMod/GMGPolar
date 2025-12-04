#pragma once

#include <cmath>

#include "../exactSolution.h"

class CartesianR2_ShafranovGeometry : public ExactSolution
{
public:
    CartesianR2_ShafranovGeometry() = default;
    explicit CartesianR2_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);
    virtual ~CartesianR2_ShafranovGeometry() = default;

    double exact_solution(double r, double theta) const override;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
