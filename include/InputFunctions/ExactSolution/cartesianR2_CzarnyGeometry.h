#pragma once

#include <cmath>

#include "../exactSolution.h"

class CartesianR2_CzarnyGeometry : public ExactSolution
{
public:
    explicit CartesianR2_CzarnyGeometry();
    explicit CartesianR2_CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon,
                                        const double& ellipticity_e);

    virtual ~CartesianR2_CzarnyGeometry() = default;

    double exact_solution(const double& r, const double& theta, const double& sin_theta,
                          const double& cos_theta) const override;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
