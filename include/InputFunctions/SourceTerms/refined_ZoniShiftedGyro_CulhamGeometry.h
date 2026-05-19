#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

namespace gmgpolar
{

class Refined_ZoniShiftedGyro_CulhamGeometry
{
public:
    explicit Refined_ZoniShiftedGyro_CulhamGeometry(PolarGrid<Kokkos::HostSpace> const& grid, double Rmax);
    KOKKOS_DEFAULTED_FUNCTION
    Refined_ZoniShiftedGyro_CulhamGeometry(const Refined_ZoniShiftedGyro_CulhamGeometry&) = default;

    KOKKOS_FUNCTION double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid<Kokkos::HostSpace> grid_;
    const double Rmax = 1.3;
};
} // namespace gmgpolar
