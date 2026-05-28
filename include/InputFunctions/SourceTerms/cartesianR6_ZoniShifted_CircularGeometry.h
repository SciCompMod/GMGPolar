#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

namespace gmgpolar
{

class CartesianR6_ZoniShifted_CircularGeometry
{
public:
    explicit CartesianR6_ZoniShifted_CircularGeometry(PolarGrid<DefaultMemorySpace> const& grid, double Rmax);
    KOKKOS_DEFAULTED_FUNCTION
    CartesianR6_ZoniShifted_CircularGeometry(const CartesianR6_ZoniShifted_CircularGeometry&) = default;

    KOKKOS_FUNCTION double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid<DefaultMemorySpace> grid_;
    const double Rmax = 1.3;
};
} // namespace gmgpolar
