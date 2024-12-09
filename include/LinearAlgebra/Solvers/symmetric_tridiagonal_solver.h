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

#include "../../common/equals.h"

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
    void solveInPlace(T* sol_rhs, T* temp1, T* temp2 = nullptr);

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const SymmetricTridiagonalSolver<U>& solver);

private:
    int matrix_dimension_;
    std::unique_ptr<T[]> main_diagonal_values_;
    std::unique_ptr<T[]> sub_diagonal_values_;
    T cyclic_corner_element_ = 0.0;
    bool is_cyclic_          = true;

    bool factorized_ = false;
    T gamma_         = 0.0;

    // Solve methods
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
    /* ---------------------------------------------------------- */
    /* Based on Cholesky Decomposition: A = L * D * L^T
    *
    * This function performs Cholesky decomposition on a 
    * symmetric tridiagonal matrix, factorizing it into 
    * a lower triangular matrix (L) and a diagonal matrix (D).
    *
    * By storing the decomposition, this approach enhances 
    * efficiency for repeated solutions, as matrix factorizations 
    * need not be recalculated each time.
    * ---------------------------------------------------------- */

    if (!factorized_) {
        for (int i = 1; i < matrix_dimension_; i++) {
            assert(!equals(main_diagonal(i - 1), 0.0));
            sub_diagonal(i - 1) /= main_diagonal(i - 1);
            main_diagonal(i) -= sub_diagonal(i - 1) * sub_diagonal(i - 1) * main_diagonal(i - 1);
        }
        factorized_ = true;
    }

    for (int i = 1; i < matrix_dimension_; i++) {
        x[i] -= sub_diagonal(i - 1) * x[i - 1];
    }

    for (int i = 0; i < matrix_dimension_; i++) {
        assert(!equals(main_diagonal(i), 0.0));
        x[i] /= main_diagonal(i);
    }

    for (int i = matrix_dimension_ - 2; i >= 0; i--) {
        x[i] -= sub_diagonal(i) * x[i + 1];
    }
}

// ------------------------- //
// Cyclic Tridiagonal Solver //
// ------------------------- //

/* 
 * This algorithm implements the Tridiagonal Matrix Algorithm (TDMA) for solving 
 * symmetric tridiagonal systems of equations, specifically designed to handle 
 * cyclic boundary conditions. The implementation is based on principles outlined 
 * in the following resource: 
 * https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm.
 */

template <typename T>
void SymmetricTridiagonalSolver<T>::solveSymmetricCyclicTridiagonal(T* x, T* u, T* scratch)
{
    /* ---------------------------------------------------------- */
    /* Cholesky Decomposition: A = L * D * L^T 
     * This step factorizes the tridiagonal matrix into lower 
     * triangular (L) and diagonal (D) matrices. While this 
     * approach may be slightly less stable, it can offer improved 
     * performance for repeated solves due to the factorization 
     * being stored internally.
     * ---------------------------------------------------------- */

    if (!factorized_) {
        gamma_ = -main_diagonal(0);
        main_diagonal(0) -= gamma_;
        main_diagonal(matrix_dimension_ - 1) -= cyclic_corner_element() * cyclic_corner_element() / gamma_;

        for (int i = 1; i < matrix_dimension_; i++) {
            sub_diagonal(i - 1) /= main_diagonal(i - 1);
            main_diagonal(i) -= sub_diagonal(i - 1) * sub_diagonal(i - 1) * main_diagonal(i - 1);
        }
        factorized_ = true;
    }

    u[0] = gamma_;
    for (int i = 1; i < matrix_dimension_; i++) {
        x[i] -= sub_diagonal(i - 1) * x[i - 1];
        if (i < matrix_dimension_ - 1)
            u[i] = 0.0 - sub_diagonal(i - 1) * u[i - 1];
        else
            u[i] = cyclic_corner_element() - sub_diagonal(i - 1) * u[i - 1];
    }

    for (int i = 0; i < matrix_dimension_; i++) {
        x[i] /= main_diagonal(i);
        u[i] /= main_diagonal(i);
    }

    for (int i = matrix_dimension_ - 2; i >= 0; i--) {
        x[i] -= sub_diagonal(i) * x[i + 1];
        u[i] -= sub_diagonal(i) * u[i + 1];
    }

    const double dot_product_x_v = x[0] + cyclic_corner_element() / gamma_ * x[matrix_dimension_ - 1];
    const double dot_product_u_v = u[0] + cyclic_corner_element() / gamma_ * u[matrix_dimension_ - 1];
    const double factor          = dot_product_x_v / (1.0 + dot_product_u_v);

    // std::cout<<u[matrix_dimension_ - 1]<<std::endl;

    for (int i = 0; i < matrix_dimension_; i++) {
        x[i] -= factor * u[i];
    }
}