#pragma once

#include <Kokkos_Core.hpp>

namespace gmgpolar
{

template <typename T, class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space>
using AllocatableVector = Kokkos::View<T*, Kokkos::LayoutRight, MemorySpace>;
template <typename T, class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space>
using Vector = Kokkos::View<T*, Kokkos::LayoutRight, MemorySpace> const;
template <typename T, class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space>
using ConstVector = Kokkos::View<const T*, Kokkos::LayoutRight, MemorySpace> const;

template <typename T>
using HostAllocatableVector = AllocatableVector<T, Kokkos::HostSpace>;
template <typename T>
using HostVector = Vector<T, Kokkos::HostSpace>;
template <typename T>
using HostConstVector = ConstVector<T, Kokkos::HostSpace>;

} // namespace gmgpolar
