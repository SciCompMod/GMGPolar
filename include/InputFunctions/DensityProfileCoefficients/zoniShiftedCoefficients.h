#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class ZoniShiftedCoefficients : public DensityProfileCoefficients
{
public:
    ZoniShiftedCoefficients() = default;
    explicit ZoniShiftedCoefficients(double Rmax, double alpha);
    virtual ~ZoniShiftedCoefficients() = default;

    double alpha(double r, double theta) const override;
    double beta(double r, double theta) const override;

    double getAlphaJump() const override;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.678 * 1.3;
};
