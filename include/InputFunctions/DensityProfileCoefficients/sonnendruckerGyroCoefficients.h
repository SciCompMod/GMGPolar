#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class SonnendruckerGyroCoefficients : public DensityProfileCoefficients { 
public:
    SonnendruckerGyroCoefficients() = default;
    explicit SonnendruckerGyroCoefficients(const double& Rmax, const double& alpha);
    virtual ~SonnendruckerGyroCoefficients() = default;

    double alpha(const double& r) const override;
    double beta(const double& r) const override;

    double getAlphaJump() const override;

private:
    const double Rmax = 1.3;
    const double alpha_jump = 0.66;
};

#include "sonnendruckerGyroCoefficients.inl"