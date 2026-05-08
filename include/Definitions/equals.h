#pragma once

#include <limits>
#include <cmath>
#include <algorithm>

#include <Kokkos_Core.hpp>

namespace gmgpolar
{

template <typename T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_floating_point_v<T>, bool> equals(T lhs, T rhs)
{
    return std::abs(lhs - rhs) <=
           (1e3 * std::numeric_limits<T>::epsilon()) * std::max(1.0, std::max(std::abs(lhs), std::abs(rhs)));
}

} // namespace gmgpolar
