#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <limits>
#include <memory>
#include <omp.h>
#include <optional>
#include <sstream>
#include <tuple>
#include <unistd.h>
#include <vector>

#include <fstream>
#include <iostream>

#include "vector.h"

template <typename T>
class DiagonalSolver
{
public:
    DiagonalSolver();
    DiagonalSolver(const DiagonalSolver& other);
    DiagonalSolver(DiagonalSolver&& other) noexcept;

    explicit DiagonalSolver(const int matrix_dimension);

    DiagonalSolver& operator=(const DiagonalSolver& other);
    DiagonalSolver& operator=(DiagonalSolver&& other) noexcept;

    int rows() const;
    int columns() const;

    const T& diagonal(const int index) const;
    T& diagonal(const int index);

    void solveInPlace(T* sol_rhs) const;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const DiagonalSolver<U>& solver);

private:
    int matrix_dimension_;
    AllocatableVector<T> diagonal_values_;
};

template <typename U>
std::ostream& operator<<(std::ostream& stream, const DiagonalSolver<U>& solver)
{
    stream << "Diagonal Matrix (Dimension: " << solver.matrix_dimension_ << " x " << solver.matrix_dimension_ << ")\n";

    stream << "Diagonal Elements: [";
    for (int i = 0; i < solver.matrix_dimension_; ++i) {
        stream << solver.diagonal(i);
        if (i != solver.matrix_dimension_ - 1)
            stream << ", ";
    }
    stream << "]\n";

    return stream;
}

// default construction
template <typename T>
DiagonalSolver<T>::DiagonalSolver()
    : matrix_dimension_(0)
{
}

// copy construction
template <typename T>
DiagonalSolver<T>::DiagonalSolver(const DiagonalSolver& other)
    : matrix_dimension_(other.matrix_dimension_)
    , diagonal_values_("Diagonal matrix values", matrix_dimension_)
{
    Kokkos::deep_copy(diagonal_values_, other.diagonal_values_);
}

// copy assignment
template <typename T>
DiagonalSolver<T>& DiagonalSolver<T>::operator=(const DiagonalSolver& other)
{
    if (this == &other) {
        // Self-assignment, no work needed
        return *this;
    }
    // Only allocate new memory if the sizes are different
    if (matrix_dimension_ != other.matrix_dimension_) {
        matrix_dimension_ = other.matrix_dimension_;
        diagonal_values_  = Vector<T>("Diagonal matrix values", matrix_dimension_);
    }
    Kokkos::deep_copy(diagonal_values_, other.diagonal_values_);
    return *this;
}

// move construction
template <typename T>
DiagonalSolver<T>::DiagonalSolver(DiagonalSolver&& other) noexcept
    : matrix_dimension_(other.matrix_dimension_)
    , diagonal_values_(std::move(other.diagonal_values_))
{
    other.matrix_dimension_ = 0;
}

// move assignment
template <typename T>
DiagonalSolver<T>& DiagonalSolver<T>::operator=(DiagonalSolver&& other) noexcept
{
    matrix_dimension_       = other.matrix_dimension_;
    diagonal_values_        = std::move(other.diagonal_values_);
    other.matrix_dimension_ = 0;
    return *this;
}

template <typename T>
DiagonalSolver<T>::DiagonalSolver(const int matrix_dimension)
    : matrix_dimension_(matrix_dimension)
    , diagonal_values_("Diagonal matrix values", matrix_dimension_)
{
    assert(matrix_dimension_ >= 1);
    Kokkos::deep_copy(diagonal_values_, T(0));
}

template <typename T>
int DiagonalSolver<T>::rows() const
{
    assert(matrix_dimension_ >= 0);
    return matrix_dimension_;
}
template <typename T>
int DiagonalSolver<T>::columns() const
{
    assert(matrix_dimension_ >= 0);
    return matrix_dimension_;
}

template <typename T>
const T& DiagonalSolver<T>::diagonal(const int index) const
{
    assert(index >= 0);
    assert(index < matrix_dimension_);
    return diagonal_values_(index);
}
template <typename T>
T& DiagonalSolver<T>::diagonal(const int index)
{
    assert(index >= 0);
    assert(index < matrix_dimension_);
    return diagonal_values_(index);
}

// --------------- //
// Diagonal Solver //
// --------------- //

template <typename T>
void DiagonalSolver<T>::solveInPlace(T* sol_rhs) const
{
    for (int i = 0; i < matrix_dimension_; i++) {
        sol_rhs[i] /= diagonal(i);
    }
}
