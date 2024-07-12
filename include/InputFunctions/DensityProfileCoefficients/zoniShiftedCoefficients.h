#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class ZoniShiftedCoefficients : public DensityProfileCoefficients { 
public:
    ZoniShiftedCoefficients() = default;
    explicit ZoniShiftedCoefficients(const double& Rmax, const double& alpha);
    virtual ~ZoniShiftedCoefficients() = default;

    double alpha(const double& r) const override;
    double beta(const double& r) const override;

    double getAlphaJump() const override;

private:
    const double Rmax = 1.3;
    const double alpha_jump = 0.7081;
};

#include "zoniShiftedCoefficients.inl"