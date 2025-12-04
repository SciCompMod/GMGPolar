#pragma once

#include <cmath>

#include "../exactSolution.h"

class Refined_CzarnyGeometry : public ExactSolution
{
public:
    explicit Refined_CzarnyGeometry();
    explicit Refined_CzarnyGeometry(double Rmax, double inverse_aspect_ratio_epsilon,
                                    double ellipticity_e);

    virtual ~Refined_CzarnyGeometry() = default;

    double exact_solution(double r, double theta) const override;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
