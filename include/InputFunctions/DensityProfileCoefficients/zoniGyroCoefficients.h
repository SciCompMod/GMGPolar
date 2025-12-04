#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class ZoniGyroCoefficients : public DensityProfileCoefficients
{
public:
    ZoniGyroCoefficients() = default;
    explicit ZoniGyroCoefficients(double Rmax, double alpha);
    virtual ~ZoniGyroCoefficients() = default;

    double alpha(double r, double theta) const override;
    double beta(double r, double theta) const override;

    double getAlphaJump() const override;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.4837 * 1.3;
};
