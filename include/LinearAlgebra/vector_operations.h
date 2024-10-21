#pragma once

#include <cmath>
#include <omp.h>
#include <stdexcept>

#include "../common/equals.h"
#include "vector.h"

// clang-format off

template <typename T>
bool equals(const Vector<T>& lhs, const Vector<T>& rhs)
{
    if (lhs.size() != rhs.size())
    {
        return false;
    }

    const std::size_t n = lhs.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!equals(lhs[i], rhs[i]))
        {
            return false;
        }
    }
    return true;
}

template <typename T>
void assign(Vector<T>& lhs, const T& value)
{
    std::size_t n = lhs.size();
    #pragma omp parallel for if (n > 100'000)
    for (std::size_t i = 0; i < n; ++i)
    {
        lhs[i] = value;
    }
}

template <typename T>
void add(Vector<T>& result, const Vector<T>& x)
{
    if (result.size() != x.size())
    {
        throw std::invalid_argument("Vectors must be of the same size.");
    }
    std::size_t n = result.size();
    #pragma omp parallel for if (n > 100'000)
    for (std::size_t i = 0; i < n; ++i)
    {
        result[i] += x[i];
    }
}

template <typename T>
void subtract(Vector<T>& result, const Vector<T>& x)
{
    if (result.size() != x.size())
    {
        throw std::invalid_argument("Vectors must be of the same size.");
    }
    std::size_t n = result.size();
    #pragma omp parallel for if (n > 100'000)
    for (std::size_t i = 0; i < n; ++i)
    {
        result[i] -= x[i];
    }
}

template <typename T>
void linear_combination(Vector<T>& x, const T& alpha, const Vector<T>& y, const T& beta)
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("Vectors must be of the same size.");
    }
    std::size_t n = x.size();
    #pragma omp parallel for if (n > 100'000)
    for (std::size_t i = 0; i < n; ++i)
    {
        x[i] = alpha * x[i] + beta * y[i];
    }
}

template <typename T>
void multiply(Vector<T>& x, const T& alpha)
{
    std::size_t n = x.size();
    #pragma omp parallel for if (n > 100'000)
    for (std::size_t i = 0; i < n; ++i)
    {
        x[i] *= alpha;
    }
}

template <typename T>
T dot_product(const Vector<T>& lhs, const Vector<T>& rhs)
{
    if (lhs.size() != rhs.size())
    {
        throw std::invalid_argument("Vectors must be of the same size.");
    }

    T result = 0.0;
    std::size_t n = lhs.size();
    #pragma omp parallel for reduction(+ : result) if (n > 100'000)
    for (std::size_t i = 0; i < n; ++i)
    {
        result += lhs[i] * rhs[i];
    }
    return result;
}

template <typename T>
T l1_norm(const Vector<T>& x)
{
    T result = 0.0;
    std::size_t n = x.size();
    #pragma omp parallel for reduction(+ : result) if (n > 100'000)
    for (std::size_t i = 0; i < n; ++i)
    {
        result += std::abs(x[i]);
    }
    return result;
}

template <typename T>
T l2_norm_squared(const Vector<T>& x)
{
    T result = 0.0;
    std::size_t n = x.size();
    #pragma omp parallel for reduction(+ : result) if (n > 100'000)
    for (std::size_t i = 0; i < n; ++i)
    {
        result += x[i] * x[i];
    }
    return result;
}

template <typename T>
T infinity_norm(const Vector<T>& x)
{
    T result = 0.0;
    std::size_t n = x.size();
    #pragma omp parallel for reduction(max : result) if (n > 100'000)
    for (std::size_t i = 0; i < n; ++i)
    {
        T abs_value = std::abs(x[i]);
        if (abs_value > result)
        {
            result = abs_value;
        }
    }
    return result;
}

// clang-format on