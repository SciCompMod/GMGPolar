#pragma once

namespace concepts
{

template <typename T>
concept BoundaryConditions = requires(const T bcs, double r, double theta) {
    { bcs.u_D(r, theta) } -> std::convertible_to<double>;
    { bcs.u_D_Interior(r, theta) } -> std::convertible_to<double>;
};

} // namespace concepts
