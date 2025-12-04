#pragma once

#include <cmath>

#include "../densityProfileCoefficients.h"

class SonnendruckerCoefficients : public DensityProfileCoefficients
{
public:
    SonnendruckerCoefficients() = default;
    explicit SonnendruckerCoefficients(double Rmax, double alpha);
    virtual ~SonnendruckerCoefficients() = default;

    double alpha(double r, double theta) const override;
    double beta(double r, double theta) const override;

    double getAlphaJump() const override;

private:
    const double Rmax       = 1.3;
    const double alpha_jump = 0.66 * 1.3;
};
