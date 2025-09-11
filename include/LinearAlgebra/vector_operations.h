#pragma once

#include <cmath>
#include <omp.h>
#include <stdexcept>

#include "../common/equals.h"
#include "vector.h"
#include <Kokkos_Core.hpp>

template <typename T>
bool equals(const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> lhs,
            const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs)
{
    if (lhs.extent(0) != rhs.extent(0)) {
        return false;
    }

    const std::size_t n = lhs.extent(0);
    for (std::size_t i = 0; i < n; ++i) {
        if (!equals(lhs[i], rhs[i])) {
            return false;
        }
    }
    return true;
}

template <typename T>
void assign(Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> lhs, const T& value)
{
    std::size_t n = lhs.extent(0);
#pragma omp parallel for if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        lhs[i] = value;
    }
}

template <typename T>
void add(Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> result,
         const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> x)
{
    if (result.extent(0) != x.extent(0)) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }
    std::size_t n = result.extent(0);
#pragma omp parallel for if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        result[i] += x(i);
    }
}

template <typename T>
void subtract(Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> result,
              const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> x)
{
    if (result.extent(0) != x.extent(0)) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }
    std::size_t n = result.extent(0);
#pragma omp parallel for if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        result[i] -= x(i);
    }
}

template <typename T>
void linear_combination(Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> x, const T& alpha,
                        const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> y, const T& beta)
{
    if (x.extent(0) != y.extent(0)) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }
    std::size_t n = x.extent(0);
#pragma omp parallel for if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        x(i) = alpha * x(i) + beta * y[i];
    }
}

template <typename T>
void multiply(Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> x, const T& alpha)
{
    std::size_t n = x.extent(0);
#pragma omp parallel for if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        x(i) *= alpha;
    }
}

template <typename T>
T dot_product(const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> lhs,
              const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs)
{
    if (lhs.extent(0) != rhs.extent(0)) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }

    T result      = 0.0;
    std::size_t n = lhs.extent(0);
#pragma omp parallel for reduction(+ : result) if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        result += lhs[i] * rhs[i];
    }
    return result;
}

template <typename T>
T l1_norm(const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> x)
{
    T result      = 0.0;
    std::size_t n = x.extent(0);
#pragma omp parallel for reduction(+ : result) if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        result += std::abs(x(i));
    }
    return result;
}

template <typename T>
T l2_norm_squared(const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> x)
{
    T result      = 0.0;
    std::size_t n = x.extent(0);
#pragma omp parallel for reduction(+ : result) if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        result += x(i) * x(i);
    }
    return result;
}

template <typename T>
T l2_norm(const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> x)
{
    const std::size_t n = x.extent(0);
    // 1) find the largest absolute value
    T scale = 0.0;
#pragma omp parallel for reduction(max : scale) if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        T abs_val = std::abs(x(i));
        if (abs_val > scale) {
            scale = abs_val;
        }
    }
    // 2) accumulate sum of squares of scaled entries
    T sum = 0.0;
#pragma omp parallel for reduction(+ : sum) if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        T value = x(i) / scale;
        sum +=value * value;
    }
    // 3) rescale
    return  scale *std::sqrt(sum);
}

template <typename T>
T infinity_norm(const Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace> x)
{
    T result      = 0.0;
    std::size_t n = x.extent(0);
#pragma omp parallel for reduction(max : result) if (n > 10'000)
    for (std::size_t i = 0; i < n; ++i) {
        T abs_value = std::abs(x(i));
        if (abs_value > result) {
            result = abs_value;
        }
    }
    return result;
}
