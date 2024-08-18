#pragma once

class DensityProfileCoefficients {
public:
    DensityProfileCoefficients() = default;
    virtual ~DensityProfileCoefficients() = default;

    virtual double alpha(const double& r) const = 0;
    virtual double beta(const double& r) const = 0;

    virtual double getAlphaJump() const = 0;
};