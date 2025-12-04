#pragma once

#include <cmath>

#include "../exactSolution.h"

class PolarR6_ShafranovGeometry : public ExactSolution
{
public:
    PolarR6_ShafranovGeometry() = default;
    explicit PolarR6_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);
    virtual ~PolarR6_ShafranovGeometry() = default;

    double exact_solution(double r, double theta) const override;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
