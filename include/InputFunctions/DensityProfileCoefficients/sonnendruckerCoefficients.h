#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class SonnendruckerCoefficients : public DensityProfileCoefficients
{
public:
    SonnendruckerCoefficients() = default;
    explicit SonnendruckerCoefficients(const double& Rmax, const double& alpha);
    virtual ~SonnendruckerCoefficients() = default;

    double alpha(const double& r) const override;
    double beta(const double& r) const override;

    double getAlphaJump() const override;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.66 * 1.3;
};
