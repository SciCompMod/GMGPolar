#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class PoissonCoefficients : public DensityProfileCoefficients
{
public:
    PoissonCoefficients() = default;
    explicit PoissonCoefficients(double Rmax, double alpha);
    virtual ~PoissonCoefficients() = default;

    double alpha(double r, double theta) const override;
    double beta(double r, double theta) const override;

    double getAlphaJump() const override;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.5 * 1.3;
};
