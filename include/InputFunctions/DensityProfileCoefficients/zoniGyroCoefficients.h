#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class ZoniGyroCoefficients : public DensityProfileCoefficients { 
public:
    ZoniGyroCoefficients() = default;
    explicit ZoniGyroCoefficients(const double& Rmax, const double& alpha);
    virtual ~ZoniGyroCoefficients() = default;

    double alpha(const double& r) const override;
    double beta(const double& r) const override;

    double getAlphaJump() const override;

private:
    const double Rmax = 1.3;
    const double alpha_jump = 0.4837;
};

#include "zoniGyroCoefficients.inl"