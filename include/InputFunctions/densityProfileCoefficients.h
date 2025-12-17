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

namespace concepts
{

template <typename T>
concept DensityProfileCoefficients =
    !std::same_as<T, DensityProfileCoefficients> && requires(const T coeffs, double r, double theta) {
        { coeffs.alpha(r, theta) } -> std::convertible_to<double>;
        { coeffs.beta(r, theta) } -> std::convertible_to<double>;
        { coeffs.getAlphaJump() } -> std::convertible_to<double>;
    };

} // namespace concepts
