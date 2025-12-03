#pragma once

#include <cmath>

#include "../exactSolution.h"

class PolarR6_CzarnyGeometry : public ExactSolution
{
public:
    explicit PolarR6_CzarnyGeometry();
    explicit PolarR6_CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon,
                                    const double& ellipticity_e);

    virtual ~PolarR6_CzarnyGeometry() = default;

    double exact_solution(const double& r, const double& theta) const override;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
