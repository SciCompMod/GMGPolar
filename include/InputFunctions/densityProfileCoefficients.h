#pragma once

class DensityProfileCoefficients
{
public:
    DensityProfileCoefficients()          = default;
    virtual ~DensityProfileCoefficients() = default;

    virtual double alpha(double r, double theta) const = 0;
    virtual double beta(double r, double theta) const  = 0;

    // Only used in custom mesh generation -> refinement_radius
    virtual double getAlphaJump() const = 0;
};