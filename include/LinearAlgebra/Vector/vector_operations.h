#pragma once

#include <cmath>
#include <omp.h>
#include <stdexcept>

#include <Kokkos_Core.hpp>

#include "vector.h"
#include "../../Definitions/equals.h"

namespace gmgpolar
{

template <typename T>
bool equals(ConstVector<T> lhs, ConstVector<T> rhs)
{
    if (lhs.size() != rhs.size()) {
        return false;
    }

    const std::size_t n = lhs.size();
    int mismatches      = 0;
    Kokkos::parallel_reduce(
        "equals", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(const std::size_t i, int& local_mismatches) {
            if (!equals(lhs(i), rhs(i))) {
                ++local_mismatches;
            }
        },
        mismatches);
    return mismatches == 0;
}

template <typename T>
void assign(Vector<T> lhs, const T& value)
{
    const std::size_t n = lhs.size();
    Kokkos::parallel_for("assign", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const std::size_t i) { lhs(i) = value; });
}

template <typename T>
void add(Vector<T> result, ConstVector<T> x)
{
    if (result.size() != x.size()) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }
    const std::size_t n = result.size();
    Kokkos::parallel_for("add", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const std::size_t i) { result(i) += x(i); });
}

template <typename T>
void subtract(Vector<T> result, ConstVector<T> x)
{
    if (result.size() != x.size()) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }
    const std::size_t n = result.size();
    Kokkos::parallel_for(
        "subtract", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const std::size_t i) { result(i) -= x(i); });
}

template <typename T>
void linear_combination(Vector<T> x, const T& alpha, ConstVector<T> y, const T& beta)
{
    if (x.size() != y.size()) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }
    const std::size_t n = x.size();
    Kokkos::parallel_for(
        "linear_combination", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(const std::size_t i) { x(i) = alpha * x(i) + beta * y(i); });
}

template <typename T>
void multiply(Vector<T> x, const T& alpha)
{
    const std::size_t n = x.size();
    Kokkos::parallel_for(
        "multiply", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const std::size_t i) { x(i) *= alpha; });
}

template <typename T>
T dot_product(ConstVector<T> lhs, ConstVector<T> rhs)
{
    if (lhs.size() != rhs.size()) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }
    const std::size_t n = lhs.size();
    T result            = T{0};
    Kokkos::parallel_reduce(
        "dot_product", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(const std::size_t i, T& local_sum) { local_sum += lhs(i) * rhs(i); }, result);
    return result;
}

template <typename T>
T l1_norm(ConstVector<T> x)
{
    const std::size_t n = x.size();
    T result            = T{0};
    Kokkos::parallel_reduce(
        "l1_norm", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(const std::size_t i, T& local_sum) { local_sum += Kokkos::abs(x(i)); }, result);
    return result;
}

// Underflow- and overflow-resistant implementation of the L2 norm
template <typename T>
T l2_norm(ConstVector<T> x)
{
    const std::size_t n = x.size();

    // 1) Find the largest absolute value
    T scale = T{0};
    Kokkos::parallel_reduce(
        "l2_norm_scale", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(const std::size_t i, T& local_max) {
            const T abs_val = Kokkos::abs(x(i));
            if (abs_val > local_max) {
                local_max = abs_val;
            }
        },
        Kokkos::Max<T>(scale));

    // 2) If the largest absolute value is zero, the norm is zero
    if (scale == T{0})
        return T{0};

    // 3) Accumulate sum of squares of scaled entries
    T sum = T{0};
    Kokkos::parallel_reduce(
        "l2_norm_sum", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(const std::size_t i, T& local_sum) {
            const T value = x(i) / scale;
            local_sum += value * value;
        },
        sum);

    // 4) Rescale
    return scale * Kokkos::sqrt(sum);
}

template <typename T>
T infinity_norm(ConstVector<T> x)
{
    const std::size_t n = x.size();
    T result            = T{0};
    Kokkos::parallel_reduce(
        "infinity_norm", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(const std::size_t i, T& local_max) {
            const T abs_value = Kokkos::abs(x(i));
            if (abs_value > local_max) {
                local_max = abs_value;
            }
        },
        Kokkos::Max<T>(result));
    return result;
}

} // namespace gmgpolar
