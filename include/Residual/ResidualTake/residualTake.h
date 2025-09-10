#pragma once

#include "../residual.h"

class ResidualTake : public Residual
{
public:
    explicit ResidualTake(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                          const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior,
                          const int num_omp_threads);
    ~ResidualTake() override = default;

    void computeResidual(Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> result,
                         const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs,
                         const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x) const override;

private:
    void applyCircleSection(const int i_r, Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> result,
                            const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs,
                            const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x) const;
    void applyRadialSection(const int i_theta, Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> result,
                            const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs,
                            const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x) const;
};