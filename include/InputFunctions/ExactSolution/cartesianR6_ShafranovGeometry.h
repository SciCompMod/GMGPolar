#pragma once

#include <cmath>

#include "../exactSolution.h"

class CartesianR6_ShafranovGeometry : public ExactSolution { 
public:
    CartesianR6_ShafranovGeometry() = default;
    explicit CartesianR6_ShafranovGeometry(const double& Rmax, const double& elongation_kappa, const double& shift_delta);
    virtual ~CartesianR6_ShafranovGeometry() = default;

    double exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta = 0.2;
};
