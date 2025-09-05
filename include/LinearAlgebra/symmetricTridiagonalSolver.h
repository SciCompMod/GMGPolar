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
 * Solver for symmetric tridiagonal systems,
 * with support for the cyclic (periodic) case.
 * 
* The algorithm uses an in-place Cholesky factorization specialized
 * for tridiagonal matrices. Once factorized, the same matrix can
 * be solved efficiently for multiple right-hand sides without
 * recomputing the decomposition.
 *
 * Usage:
 * - Construct with the system dimension (and optionally enable cyclic mode).
 * - Set the main diagonal, sub-diagonal, and, if cyclic, the corner element.
 * - Call the appropriate solve method with the right-hand side vector.
 * - Provide a scratch 'temp' vector of size matrix_dimension_ in the cyclic case.
 */
template <typename T>
class SymmetricTridiagonalSolver
{
public:
    SymmetricTridiagonalSolver();
    SymmetricTridiagonalSolver(const SymmetricTridiagonalSolver& other);
    SymmetricTridiagonalSolver(SymmetricTridiagonalSolver&& other) noexcept;

    explicit SymmetricTridiagonalSolver(const int matrix_dimension);

    SymmetricTridiagonalSolver& operator=(const SymmetricTridiagonalSolver& other);
    SymmetricTridiagonalSolver& operator=(SymmetricTridiagonalSolver&& other) noexcept;

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
    void solveInPlace(T* sol_rhs, T* temp = nullptr);

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const SymmetricTridiagonalSolver<U>& solver);

private:
    int matrix_dimension_;
    std::unique_ptr<T[]> main_diagonal_values_;
    std::unique_ptr<T[]> sub_diagonal_values_;
    T cyclic_corner_element_ = 0.0;
    bool is_cyclic_          = true;

    bool factorized_ = false;
    T gamma_         = 0.0; // Used in Shermann-Morrison factorization A = B + u*v^T

    // Solve methods
    void solveSymmetricTridiagonal(T* x);
    void solveSymmetricCyclicTridiagonal(T* x, T* u);
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
{
}

// copy construction
template <typename T>
SymmetricTridiagonalSolver<T>::SymmetricTridiagonalSolver(const SymmetricTridiagonalSolver& other)
    : matrix_dimension_(other.matrix_dimension_)
    , main_diagonal_values_(std::make_unique<T[]>(matrix_dimension_))
    , sub_diagonal_values_(std::make_unique<T[]>(matrix_dimension_ - 1))
    , cyclic_corner_element_(other.cyclic_corner_element_)
    , is_cyclic_(other.is_cyclic_)
{
    std::copy(other.main_diagonal_values_.get(), other.main_diagonal_values_.get() + matrix_dimension_,
              main_diagonal_values_.get());
    std::copy(other.sub_diagonal_values_.get(), other.sub_diagonal_values_.get() + matrix_dimension_ - 1,
              sub_diagonal_values_.get());
}

// copy assignment
template <typename T>
SymmetricTridiagonalSolver<T>& SymmetricTridiagonalSolver<T>::operator=(const SymmetricTridiagonalSolver& other)
{
    if (this == &other) {
        // Self-assignment, no work needed
        return *this;
    }
    // Only allocate new memory if the sizes are different
    if (matrix_dimension_ != other.matrix_dimension_) {
        matrix_dimension_     = other.matrix_dimension_;
        main_diagonal_values_ = std::make_unique<T[]>(matrix_dimension_);
        sub_diagonal_values_  = std::make_unique<T[]>(matrix_dimension_ - 1);
    }
    cyclic_corner_element_ = other.cyclic_corner_element_;
    is_cyclic_             = other.is_cyclic_;
    std::copy(other.main_diagonal_values_.get(), other.main_diagonal_values_.get() + matrix_dimension_,
              main_diagonal_values_.get());
    std::copy(other.sub_diagonal_values_.get(), other.sub_diagonal_values_.get() + matrix_dimension_ - 1,
              sub_diagonal_values_.get());
    return *this;
}

// move construction
template <typename T>
SymmetricTridiagonalSolver<T>::SymmetricTridiagonalSolver(SymmetricTridiagonalSolver&& other) noexcept
    : matrix_dimension_(other.matrix_dimension_)
    , main_diagonal_values_(std::move(other.main_diagonal_values_))
    , sub_diagonal_values_(std::move(other.sub_diagonal_values_))
    , cyclic_corner_element_(other.cyclic_corner_element_)
    , is_cyclic_(other.is_cyclic_)
{
    other.matrix_dimension_      = 0;
    other.cyclic_corner_element_ = 0.0;
    other.is_cyclic_             = true;
}

// move assignment
template <typename T>
SymmetricTridiagonalSolver<T>& SymmetricTridiagonalSolver<T>::operator=(SymmetricTridiagonalSolver&& other) noexcept
{
    matrix_dimension_            = other.matrix_dimension_;
    main_diagonal_values_        = std::move(other.main_diagonal_values_);
    sub_diagonal_values_         = std::move(other.sub_diagonal_values_);
    cyclic_corner_element_       = other.cyclic_corner_element_;
    is_cyclic_                   = other.is_cyclic_;
    other.matrix_dimension_      = 0;
    other.cyclic_corner_element_ = 0.0;
    other.is_cyclic_             = true;
    return *this;
}

template <typename T>
SymmetricTridiagonalSolver<T>::SymmetricTridiagonalSolver(const int matrix_dimension)
    : matrix_dimension_(matrix_dimension)
    , main_diagonal_values_(std::make_unique<T[]>(matrix_dimension_))
    , sub_diagonal_values_(std::make_unique<T[]>(matrix_dimension_ - 1))
    , cyclic_corner_element_(0.0)
    , is_cyclic_(true)
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
    assert(this->matrix_dimension_ >= 0);
    return this->matrix_dimension_;
}
template <typename T>
int SymmetricTridiagonalSolver<T>::columns() const
{
    assert(this->matrix_dimension_ >= 0);
    return this->matrix_dimension_;
}

template <typename T>
const T& SymmetricTridiagonalSolver<T>::main_diagonal(const int index) const
{
    assert(index >= 0);
    assert(index < this->matrix_dimension_);
    return this->main_diagonal_values_[index];
}
template <typename T>
T& SymmetricTridiagonalSolver<T>::main_diagonal(const int index)
{
    assert(index >= 0);
    assert(index < this->matrix_dimension_);
    return this->main_diagonal_values_[index];
}

template <typename T>
const T& SymmetricTridiagonalSolver<T>::sub_diagonal(const int index) const
{
    assert(index >= 0);
    assert(index < this->matrix_dimension_ - 1);
    return this->sub_diagonal_values_[index];
}
template <typename T>
T& SymmetricTridiagonalSolver<T>::sub_diagonal(const int index)
{
    assert(index >= 0);
    assert(index < this->matrix_dimension_ - 1);
    return this->sub_diagonal_values_[index];
}

template <typename T>
const T& SymmetricTridiagonalSolver<T>::cyclic_corner_element() const
{
    assert(is_cyclic_);
    return this->cyclic_corner_element_;
}
template <typename T>
T& SymmetricTridiagonalSolver<T>::cyclic_corner_element()
{
    assert(is_cyclic_);
    return this->cyclic_corner_element_;
}

template <typename T>
void SymmetricTridiagonalSolver<T>::solveInPlace(T* sol_rhs, T* temp)
{
    assert(matrix_dimension_ >= 2);
    assert(sol_rhs != nullptr);
    if (is_cyclic_) {
        assert(temp != nullptr);
        solveSymmetricCyclicTridiagonal(sol_rhs, temp);
    }
    else {
        solveSymmetricTridiagonal(sol_rhs);
    }
}

/* ----------------------------------------------------------
 * This algorithm implements the Tridiagonal Matrix Algorithm (TDMA) for solving 
 * symmetric tridiagonal systems of equations. The implementation is based on 
 * the principles outlined in the following resource: 
 * https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm.
 * 
 * It is based on the Cholesky Decomposition: A = L * D * L^T
 *
 * This function performs Cholesky decomposition on a 
 * symmetric tridiagonal matrix, factorizing it into 
 * a lower triangular matrix (L) and a diagonal matrix (D).
 *
 * By storing the decomposition, this approach enhances 
 * efficiency for repeated solutions, as matrix factorizations 
 * need not be recalculated each time.
 */
template <typename T>
void SymmetricTridiagonalSolver<T>::solveSymmetricTridiagonal(T* x)
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

/* ----------------------------------------------------------
 * This algorithm solves cyclic symmetric tridiagonal systems.
 *
 * The system Ax = b, where A is cyclic, can be written as
 *    A = B + u v^T
 * with B tridiagonal. By the Sherman–Morrison formula:
 *    A⁻¹ b = B⁻¹ b − (B⁻¹ u) (v^T B⁻¹ b) / (1 + v^T B⁻¹ u).
 *
 * Thus, we only need to solve two standard tridiagonal systems:
 *    B y = b   and   B z = u.
 * The final solution is then
 *    x = y − z (v^T y) / (1 + v^T z).
 *
 * This reduces the cyclic case to two tridiagonal solves,
 * followed by a scalar correction step.
 */
template <typename T>
void SymmetricTridiagonalSolver<T>::solveSymmetricCyclicTridiagonal(T* x, T* u)
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
