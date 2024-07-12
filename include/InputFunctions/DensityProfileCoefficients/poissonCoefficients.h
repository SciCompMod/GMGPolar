#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class PoissonCoefficients : public DensityProfileCoefficients { 
public:
    PoissonCoefficients() = default;
    explicit PoissonCoefficients(const double& Rmax, const double& alpha);
    virtual ~PoissonCoefficients() = default;

    double alpha(const double& r) const override;
    double beta(const double& r) const override;

    double getAlphaJump() const override;

private:
    const double Rmax = 1.3;
    const double alpha_jump = 0.5;
};

#include "poissonCoefficients.inl"