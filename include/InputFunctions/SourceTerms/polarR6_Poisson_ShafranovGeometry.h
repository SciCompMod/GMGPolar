#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

namespace gmgpolar
{

class PolarR6_Poisson_ShafranovGeometry
{
public:
    explicit PolarR6_Poisson_ShafranovGeometry(PolarGrid<Kokkos::HostSpace> const& grid, double Rmax, double elongation_kappa,
                                               double shift_delta);
    KOKKOS_DEFAULTED_FUNCTION PolarR6_Poisson_ShafranovGeometry(const PolarR6_Poisson_ShafranovGeometry&) = default;

    KOKKOS_FUNCTION double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid<Kokkos::HostSpace> grid_;
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
} // namespace gmgpolar
