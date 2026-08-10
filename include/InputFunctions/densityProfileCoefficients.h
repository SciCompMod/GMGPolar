#pragma once

namespace concepts
{

template <typename T>
concept DensityProfileCoefficients = requires(const T coeffs, int i_r, int i_theta) {
    { coeffs.alpha(i_r, i_theta) } -> std::convertible_to<double>;
    { coeffs.beta(i_r, i_theta) } -> std::convertible_to<double>;
    { coeffs.getAlphaJump() } -> std::convertible_to<double>;
};

} // namespace concepts
