#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class ZoniShiftedGyroCoefficients : public DensityProfileCoefficients
{
public:
    ZoniShiftedGyroCoefficients() = default;
    explicit ZoniShiftedGyroCoefficients(double Rmax, double alpha);
    virtual ~ZoniShiftedGyroCoefficients() = default;

    double alpha(double r, double theta) const override;
    double beta(double r, double theta) const override;

    double getAlphaJump() const override;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.678 * 1.3;
};
