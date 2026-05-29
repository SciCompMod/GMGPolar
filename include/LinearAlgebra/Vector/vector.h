#pragma once

#include <Kokkos_Core.hpp>

namespace gmgpolar
{

using DefaultMemorySpace = Kokkos::DefaultExecutionSpace::memory_space;

template <typename T, class MemorySpace = DefaultMemorySpace>
using AllocatableVector = Kokkos::View<T*, Kokkos::LayoutRight, MemorySpace>;
template <typename T, class MemorySpace = DefaultMemorySpace>
using Vector = Kokkos::View<T*, Kokkos::LayoutRight, MemorySpace> const;
template <typename T, class MemorySpace = DefaultMemorySpace>
using ConstVector = Kokkos::View<const T*, Kokkos::LayoutRight, MemorySpace> const;

template <typename T>
using HostAllocatableVector = AllocatableVector<T, Kokkos::HostSpace>;
template <typename T>
using HostVector = Vector<T, Kokkos::HostSpace>;
template <typename T>
using HostConstVector = ConstVector<T, Kokkos::HostSpace>;

} // namespace gmgpolar
