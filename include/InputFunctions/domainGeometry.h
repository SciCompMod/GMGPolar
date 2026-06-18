#pragma once

/**
 * @class DomainGeometry
 * @brief A concept representing the geometric properties of a domain.
 *
 * These classes provide an interface for calculating geometric transformations and their derivatives
 * for a domain in polar coordinates (r, θ). They include functions to compute the
 * Cartesian coordinates (Fx, Fy) from polar coordinates, as well as their partial derivatives with
 * respect to r and θ.
 */

namespace concepts
{

template <typename T>
concept DomainGeometry = requires(const T geometry, double r, double theta) {
    { geometry.Fx(r, theta) } -> std::convertible_to<double>;
    { geometry.Fy(r, theta) } -> std::convertible_to<double>;
    { geometry.dFx_dr(r, theta) } -> std::convertible_to<double>;
    { geometry.dFy_dr(r, theta) } -> std::convertible_to<double>;
    { geometry.dFx_dtheta(r, theta) } -> std::convertible_to<double>;
    { geometry.dFy_dtheta(r, theta) } -> std::convertible_to<double>;
};

} // namespace concepts
