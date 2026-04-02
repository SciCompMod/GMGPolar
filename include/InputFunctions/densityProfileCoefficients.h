#pragma once

namespace concepts
{

template <typename T>
concept DensityProfileCoefficients =
    requires(const T coeffs, double r, double theta) {
        { coeffs.alpha(r, theta) } -> std::convertible_to<double>;
        { coeffs.beta(r, theta) } -> std::convertible_to<double>;
        { coeffs.getAlphaJump() } -> std::convertible_to<double>;
    };

} // namespace concepts
