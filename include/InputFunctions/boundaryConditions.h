#pragma once

class BoundaryConditions
{
public:
    BoundaryConditions()          = default;
    virtual ~BoundaryConditions() = default;

    virtual double u_D(double r, double theta) const = 0;
    // Only used if DirBC_Interior = true
    virtual double u_D_Interior(double r, double theta) const = 0;
};

namespace concepts
{

template <typename T>
concept BoundaryConditions = !std::same_as<T, BoundaryConditions> && requires(const T bcs, double r, double theta) {
    { bcs.u_D(r, theta) } -> std::convertible_to<double>;
    { bcs.u_D_Interior(r, theta) } -> std::convertible_to<double>;
};

} // namespace concepts
