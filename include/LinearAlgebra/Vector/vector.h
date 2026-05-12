#pragma once

#include <Kokkos_Core.hpp>

namespace gmgpolar
{

template <typename T, class MemorySpace = Kokkos::HostSpace>
using AllocatableVector = Kokkos::View<T*, Kokkos::LayoutRight, MemorySpace>;
template <typename T, class MemorySpace = Kokkos::HostSpace>
using Vector = Kokkos::View<T*, Kokkos::LayoutRight, MemorySpace> const;
template <typename T, class MemorySpace = Kokkos::HostSpace>
using ConstVector = Kokkos::View<const T*, Kokkos::LayoutRight, MemorySpace> const;

} // namespace gmgpolar
