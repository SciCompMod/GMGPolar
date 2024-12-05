#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class ZoniShiftedGyroCoefficients : public DensityProfileCoefficients
{
public:
    ZoniShiftedGyroCoefficients() = default;
    explicit ZoniShiftedGyroCoefficients(const double& Rmax, const double& alpha);
    virtual ~ZoniShiftedGyroCoefficients() = default;

    double alpha(const double& r) const override;
    double beta(const double& r) const override;

    double getAlphaJump() const override;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.7081 * 1.3;
};
