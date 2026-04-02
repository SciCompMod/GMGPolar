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
concept DomainGeometry = requires(const T geom, double r, double theta) {
    { geom.Fx(r, theta) } -> std::convertible_to<double>;
    { geom.Fy(r, theta) } -> std::convertible_to<double>;
    { geom.dFx_dr(r, theta) } -> std::convertible_to<double>;
    { geom.dFy_dr(r, theta) } -> std::convertible_to<double>;
    { geom.dFx_dt(r, theta) } -> std::convertible_to<double>;
    { geom.dFy_dt(r, theta) } -> std::convertible_to<double>;
};

} // namespace concepts
