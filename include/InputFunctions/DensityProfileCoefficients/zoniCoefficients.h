#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class ZoniCoefficients : public DensityProfileCoefficients
{
public:
    ZoniCoefficients() = default;
    explicit ZoniCoefficients(const double& Rmax, const double& alpha);
    virtual ~ZoniCoefficients() = default;

    double alpha(const double& r) const override;
    double beta(const double& r) const override;

    double getAlphaJump() const override;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.4837 * 1.3;
};
