#pragma once

#include <Kokkos_Core.hpp>

namespace gmgpolar
{

template <typename T>
using AllocatableVector = Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace>;
template <typename T>
using Vector = Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> const;
template <typename T>
using ConstVector = Kokkos::View<const T*, Kokkos::LayoutRight, Kokkos::HostSpace> const;

} // namespace gmgpolar
