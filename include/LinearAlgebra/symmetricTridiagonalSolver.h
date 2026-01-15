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

#include "../common/equals.h"

/*
 * SymmetricTridiagonalSolver is a class for solving symmetric tridiagonal systems of linear equations.
 * The system is represented by a matrix that has non-zero values only on the main diagonal,
 * the sub-diagonal, and (optionally) the cyclic corner element (if cyclic boundary conditions are used).
 * This class provides efficient solvers for both cyclic and non-cyclic symmetric tridiagonal matrices.
 * 
 * The class supports the following:
 * - Solving the system in-place using the `solveInPlace` method.
 * - Handling both cyclic and non-cyclic boundary conditions.
 * - Storing the matrix's main diagonal, sub-diagonal, and an optional cyclic corner element.
 * - Peforming the Cholesky-Decomposition in-place.
 * 
 * The primary method for solving the system is `solveInPlace`, which computes the solution to the system 
 * in place, updating the provided solution vector (`sol_rhs`) using intermediate storage (`temp1`, `temp2`).
 * 
 * Temporary storage (`temp1` and `temp2`) of length 'matrix_dimension_' must be provided by the user. 
 * These temporary vectors are used for intermediate calculations during the solving process:
 * - `temp1` is used in both non-cyclic and cyclic systems.
 * - `temp2` is only required for cyclic systems.
 * 
 * Usage:
 * - Instantiate the solver with a specified matrix dimension.
 * - Optionally set the cyclic boundary condition flag.
 * - Call `solveInPlace` to solve the system, passing the solution vector and the appropriate temporary storage.
 * 
 * The solver can handle both cyclic and non-cyclic matrices, and it uses efficient algorithms
 * optimized for symmetric tridiagonal systems.
 */

template <typename T>
class SymmetricTridiagonalSolver
{
public:
    SymmetricTridiagonalSolver();
    explicit SymmetricTridiagonalSolver(const int matrix_dimension);

    void is_cyclic(bool value);
    bool is_cyclic() const;

    int rows() const;
    int columns() const;

    const T& main_diagonal(const int index) const;
    T& main_diagonal(const int index);

    const T& sub_diagonal(const int index) const;
    T& sub_diagonal(const int index);

    const T& cyclic_corner_element() const;
    T& cyclic_corner_element();

    // Unified Solve method
    void solveInPlace(T* sol_rhs, T* temp1, T* temp2 = nullptr);

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const SymmetricTridiagonalSolver<U>& solver);

private:
    int matrix_dimension_;
    std::unique_ptr<T[]> main_diagonal_values_;
    std::unique_ptr<T[]> sub_diagonal_values_;
    T cyclic_corner_element_;
    bool is_cyclic_;
    bool factorized_;
    T gamma_; // Used in Shermann-Morrison factorization A = B + u*v^T

    // Solve methods:
    // The notation 'u' and 'scratch' is directly taken from the implementation
    // of the Thomas Algorithm listed on wikipedia.
    // Note that the 'scratch' vector is unused, as we moved from the
    // Thomas Algorithm to the faster Cholesky Decomposition.
    void solveSymmetricTridiagonal(T* x, T* scratch);
    void solveSymmetricCyclicTridiagonal(T* x, T* u, T* scratch);
};

template <typename U>
std::ostream& operator<<(std::ostream& stream, const SymmetricTridiagonalSolver<U>& solver)
{
    stream << "Symmetric Tridiagonal Matrix (Dimension: " << solver.matrix_dimension_ << ")\n";

    if (solver.factorized_) {
        // Print the L, D decomposition if factorized
        stream << "L Factor (Sub Diagonal Elements): [";
        for (int i = 0; i < solver.matrix_dimension_ - 1; ++i) {
            stream << solver.sub_diagonal(i);
            if (i != solver.matrix_dimension_ - 2)
                stream << ", ";
        }
        stream << "]\n";

        stream << "D Factor (Diagonal Elements): [";
        for (int i = 0; i < solver.matrix_dimension_; ++i) {
            stream << solver.main_diagonal(i);
            if (i != solver.matrix_dimension_ - 1)
                stream << ", ";
        }
        stream << "]\n";
    }
    else {
        // Print the matrix in its tridiagonal form if not factorized
        stream << "Main Diagonal: [";
        for (int i = 0; i < solver.matrix_dimension_; ++i) {
            stream << solver.main_diagonal(i);
            if (i != solver.matrix_dimension_ - 1)
                stream << ", ";
        }
        stream << "]\n";

        stream << "Sub Diagonal: [";
        for (int i = 0; i < solver.matrix_dimension_ - 1; ++i) {
            stream << solver.sub_diagonal(i);
            if (i != solver.matrix_dimension_ - 2)
                stream << ", ";
        }
        stream << "]\n";

        if (solver.is_cyclic_) {
            stream << "Cyclic Corner Element: " << solver.cyclic_corner_element() << "\n";
        }
        else {
            stream << "Matrix is not cyclic.\n";
        }
    }

    return stream;
}

// default construction
template <typename T>
SymmetricTridiagonalSolver<T>::SymmetricTridiagonalSolver()
    : matrix_dimension_(0)
    , main_diagonal_values_(nullptr)
    , sub_diagonal_values_(nullptr)
    , cyclic_corner_element_(0.0)
    , is_cyclic_(true)
    , factorized_(false)
    , gamma_(0.0)
{
}

template <typename T>
SymmetricTridiagonalSolver<T>::SymmetricTridiagonalSolver(const int matrix_dimension)
    : matrix_dimension_(matrix_dimension)
    , main_diagonal_values_(std::make_unique<T[]>(matrix_dimension_))
    , sub_diagonal_values_(std::make_unique<T[]>(matrix_dimension_ - 1))
    , cyclic_corner_element_(0.0)
    , is_cyclic_(true)
    , factorized_(false)
    , gamma_(0.0)
{
    assert(matrix_dimension_ >= 1);
    std::fill(main_diagonal_values_.get(), main_diagonal_values_.get() + matrix_dimension_, T(0));
    std::fill(sub_diagonal_values_.get(), sub_diagonal_values_.get() + matrix_dimension_ - 1, T(0));
}

template <typename T>
void SymmetricTridiagonalSolver<T>::is_cyclic(bool value)
{
    is_cyclic_ = value;
}
template <typename T>
bool SymmetricTridiagonalSolver<T>::is_cyclic() const
{
    return is_cyclic_;
}

template <typename T>
int SymmetricTridiagonalSolver<T>::rows() const
{
    return matrix_dimension_;
}
template <typename T>
int SymmetricTridiagonalSolver<T>::columns() const
{
    return matrix_dimension_;
}

template <typename T>
const T& SymmetricTridiagonalSolver<T>::main_diagonal(const int index) const
{
    assert(0 <= index && index < matrix_dimension_);
    return main_diagonal_values_[index];
}
template <typename T>
T& SymmetricTridiagonalSolver<T>::main_diagonal(const int index)
{
    assert(0 <= index && index < matrix_dimension_);
    return main_diagonal_values_[index];
}

template <typename T>
const T& SymmetricTridiagonalSolver<T>::sub_diagonal(const int index) const
{
    assert(0 <= index && index < matrix_dimension_ - 1);
    return sub_diagonal_values_[index];
}
template <typename T>
T& SymmetricTridiagonalSolver<T>::sub_diagonal(const int index)
{
    assert(0 <= index && index < matrix_dimension_ - 1);
    return sub_diagonal_values_[index];
}

template <typename T>
const T& SymmetricTridiagonalSolver<T>::cyclic_corner_element() const
{
    assert(is_cyclic_);
    return cyclic_corner_element_;
}
template <typename T>
T& SymmetricTridiagonalSolver<T>::cyclic_corner_element()
{
    assert(is_cyclic_);
    return cyclic_corner_element_;
}

template <typename T>
void SymmetricTridiagonalSolver<T>::solveInPlace(T* sol_rhs, T* temp1, T* temp2)
{
    assert(matrix_dimension_ >= 2);
    assert(sol_rhs != nullptr);
    assert(temp1 != nullptr);
    if (is_cyclic_) {
        assert(temp2 != nullptr);
        solveSymmetricCyclicTridiagonal(sol_rhs, temp1, temp2);
    }
    else {
        solveSymmetricTridiagonal(sol_rhs, temp1);
    }
}

/* 
 * This algorithm implements the Tridiagonal Matrix Algorithm (TDMA) for solving 
 * symmetric tridiagonal systems of equations. The implementation is based on 
 * the principles outlined in the following resource: 
 * https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm.
 */

template <typename T>
void SymmetricTridiagonalSolver<T>::solveSymmetricTridiagonal(T* x, T* scratch)
{
    // Cholesky Decomposition
    if (!factorized_) {
        for (int i = 1; i < matrix_dimension_; i++) {
            assert(!equals(main_diagonal(i - 1), 0.0));
            sub_diagonal(i - 1) /= main_diagonal(i - 1);
            main_diagonal(i) -= sub_diagonal(i - 1) * sub_diagonal(i - 1) * main_diagonal(i - 1);
        }
        factorized_ = true;
    }
    // Forward Substitution
    for (int i = 1; i < matrix_dimension_; i++) {
        x[i] -= sub_diagonal(i - 1) * x[i - 1];
    }
    // Diagonal Scaling
    for (int i = 0; i < matrix_dimension_; i++) {
        assert(!equals(main_diagonal(i), 0.0));
        x[i] /= main_diagonal(i);
    }
    // Backward Substitution
    for (int i = matrix_dimension_ - 2; i >= 0; i--) {
        x[i] -= sub_diagonal(i) * x[i + 1];
    }
}

// ------------------------- //
// Cyclic Tridiagonal Solver //
// ------------------------- //

template <typename T>
void SymmetricTridiagonalSolver<T>::solveSymmetricCyclicTridiagonal(T* x, T* u, T* scratch)
{
    // Cholesky Decomposition
    if (!factorized_) {
        // Shermann-Morrison Adjustment
        gamma_ = -main_diagonal(0);
        main_diagonal(0) -= gamma_;
        main_diagonal(matrix_dimension_ - 1) -= cyclic_corner_element() * cyclic_corner_element() / gamma_;

        for (int i = 1; i < matrix_dimension_; i++) {
            sub_diagonal(i - 1) /= main_diagonal(i - 1);
            main_diagonal(i) -= sub_diagonal(i - 1) * sub_diagonal(i - 1) * main_diagonal(i - 1);
        }
        factorized_ = true;
    }
    // Forward Substitution
    u[0] = gamma_;
    for (int i = 1; i < matrix_dimension_; i++) {
        x[i] -= sub_diagonal(i - 1) * x[i - 1];
        if (i < matrix_dimension_ - 1)
            u[i] = 0.0 - sub_diagonal(i - 1) * u[i - 1];
        else
            u[i] = cyclic_corner_element() - sub_diagonal(i - 1) * u[i - 1];
    }
    // Diagonal Scaling
    for (int i = 0; i < matrix_dimension_; i++) {
        x[i] /= main_diagonal(i);
        u[i] /= main_diagonal(i);
    }
    // Backward Substitution
    for (int i = matrix_dimension_ - 2; i >= 0; i--) {
        x[i] -= sub_diagonal(i) * x[i + 1];
        u[i] -= sub_diagonal(i) * u[i + 1];
    }
    // Shermann-Morrison Reonstruction
    const double dot_product_x_v = x[0] + cyclic_corner_element() / gamma_ * x[matrix_dimension_ - 1];
    const double dot_product_u_v = u[0] + cyclic_corner_element() / gamma_ * u[matrix_dimension_ - 1];
    const double factor          = dot_product_x_v / (1.0 + dot_product_u_v);

    for (int i = 0; i < matrix_dimension_; i++) {
        x[i] -= factor * u[i];
    }
}
