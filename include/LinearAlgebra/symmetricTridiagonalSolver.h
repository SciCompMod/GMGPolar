#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <tuple>
#include <vector>
#include <sstream>
#include <unistd.h>
#include <omp.h>

#include <fstream>
#include <iostream>

/* ------------------------------------------------------------------ */
/* Thomas' algorithm is not stable in general,                        */
/* but is so in several special cases, such as                        */
/* when the matrix is diagonally dominant (either by rows or columns) */
/* or symmetric positive definite.                                    */
/* ------------------------------------------------------------------ */

template<typename T>
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
    void solveInPlace(T* sol_rhs, T* temp1, T* temp2 = nullptr) const;
private:
    int matrix_dimension_;
    std::unique_ptr<T[]> main_diagonal_values_;
    std::unique_ptr<T[]> sub_diagonal_values_;
    T cyclic_corner_element_ = 0.0;
    bool is_cyclic_ = true;

    bool solver_setup_ = false;

    // Solve methods
    void solve_symmetricTridiagonal(T* x, T* scratch) const;
    void solve_symmetricCyclicTridiagonal(T* x, T* u, T* scratch) const;
};


// default construction
template<typename T>
SymmetricTridiagonalSolver<T>::SymmetricTridiagonalSolver() :
    matrix_dimension_(0),
    main_diagonal_values_(nullptr),
    sub_diagonal_values_(nullptr),
    cyclic_corner_element_(0.0),
    is_cyclic_(true)
{}

// copy construction
template<typename T>
SymmetricTridiagonalSolver<T>::SymmetricTridiagonalSolver(const SymmetricTridiagonalSolver& other) :
    matrix_dimension_(other.matrix_dimension_),
    main_diagonal_values_(std::make_unique<T[]>(matrix_dimension_)),
    sub_diagonal_values_(std::make_unique<T[]>(matrix_dimension_-1)),
    cyclic_corner_element_(other.cyclic_corner_element_),
    is_cyclic_(other.is_cyclic_)
{
    std::copy(other.main_diagonal_values_.get(), other.main_diagonal_values_.get() + matrix_dimension_, main_diagonal_values_.get());
    std::copy(other.sub_diagonal_values_.get(), other.sub_diagonal_values_.get() + matrix_dimension_-1, sub_diagonal_values_.get());
}

// copy assignment
template<typename T>
SymmetricTridiagonalSolver<T>& SymmetricTridiagonalSolver<T>::operator=(const SymmetricTridiagonalSolver& other){
    if (this == &other) {
        // Self-assignment, no work needed
        return *this;
    }
    // Only allocate new memory if the sizes are different
    if (matrix_dimension_ != other.matrix_dimension_) {
        matrix_dimension_ = other.matrix_dimension_;
        main_diagonal_values_ = std::make_unique<T[]>(matrix_dimension_);
        sub_diagonal_values_ = std::make_unique<T[]>(matrix_dimension_-1);
    }
    cyclic_corner_element_ = other.cyclic_corner_element_;
    is_cyclic_ = other.is_cyclic_;
    std::copy(other.main_diagonal_values_.get(), other.main_diagonal_values_.get() + matrix_dimension_, main_diagonal_values_.get());
    std::copy(other.sub_diagonal_values_.get(), other.sub_diagonal_values_.get() + matrix_dimension_-1, sub_diagonal_values_.get());
    return *this;
}

// move construction
template<typename T>
SymmetricTridiagonalSolver<T>::SymmetricTridiagonalSolver(SymmetricTridiagonalSolver&& other) noexcept :
    matrix_dimension_(other.matrix_dimension_),
    main_diagonal_values_(std::move(other.main_diagonal_values_)),
    sub_diagonal_values_(std::move(other.sub_diagonal_values_)),
    cyclic_corner_element_(other.cyclic_corner_element_),
    is_cyclic_(other.is_cyclic_)
{}

// move assignment
template<typename T>
SymmetricTridiagonalSolver<T>& SymmetricTridiagonalSolver<T>::operator=(SymmetricTridiagonalSolver&& other) noexcept{
    matrix_dimension_ = other.matrix_dimension_;
    main_diagonal_values_ = std::move(other.main_diagonal_values_);
    sub_diagonal_values_ = std::move(other.sub_diagonal_values_);
    cyclic_corner_element_ = other.cyclic_corner_element_;
    is_cyclic_ = other.is_cyclic_;
    return *this;
}

template<typename T>
SymmetricTridiagonalSolver<T>::SymmetricTridiagonalSolver(const int matrix_dimension):
    matrix_dimension_(matrix_dimension),
    main_diagonal_values_(std::make_unique<T[]>(matrix_dimension_)),
    sub_diagonal_values_(std::make_unique<T[]>(matrix_dimension_-1)),
    cyclic_corner_element_(0.0),
    is_cyclic_(true)
{
    assert(matrix_dimension >= 1);
}


template<typename T>
void SymmetricTridiagonalSolver<T>::is_cyclic(bool value){
    is_cyclic_ = value;
}
template<typename T>
bool SymmetricTridiagonalSolver<T>::is_cyclic() const{
    return is_cyclic_;
}

template<typename T>
int SymmetricTridiagonalSolver<T>::rows() const {
    assert(this->matrix_dimension_ >= 0);
    return this->matrix_dimension_;
}
template<typename T>
int SymmetricTridiagonalSolver<T>::columns() const {
    assert(this->matrix_dimension_ >= 0);
    return this->matrix_dimension_;
}

template<typename T>
const T& SymmetricTridiagonalSolver<T>::main_diagonal(const int index) const{
    assert(index >= 0); assert(index < this->matrix_dimension_);
    return this->main_diagonal_values_[index];
}
template<typename T>
T& SymmetricTridiagonalSolver<T>::main_diagonal(const int index){
    assert(index >= 0); assert(index < this->matrix_dimension_);
    return this->main_diagonal_values_[index];  
}

template<typename T>
const T& SymmetricTridiagonalSolver<T>::sub_diagonal(const int index) const{
    assert(index >= 0); assert(index < this->matrix_dimension_-1);
    return this->sub_diagonal_values_[index];
}
template<typename T>
T& SymmetricTridiagonalSolver<T>::sub_diagonal(const int index){
    assert(index >= 0); assert(index < this->matrix_dimension_-1);
    return this->sub_diagonal_values_[index];
}

template<typename T>
const T& SymmetricTridiagonalSolver<T>::cyclic_corner_element() const{
    assert(is_cyclic_);
    return this->cyclic_corner_element_;
}
template<typename T>
T& SymmetricTridiagonalSolver<T>::cyclic_corner_element(){
    assert(is_cyclic_);
    return this->cyclic_corner_element_;
}

template<typename T>
void SymmetricTridiagonalSolver<T>::solveInPlace(T* sol_rhs, T* temp1, T* temp2) const{
    assert(matrix_dimension_ >= 2);
    assert(sol_rhs != nullptr);
    assert(temp1 != nullptr);
    if(is_cyclic_){
        assert(temp2 != nullptr);
        solve_symmetricCyclicTridiagonal(sol_rhs, temp1, temp2);
    }
    else{
        solve_symmetricTridiagonal(sol_rhs, temp1);
    }
}

// ------------------ //
// Tridiagonal Solver //
// ------------------ //

/* Algorithm based on: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm */
template<typename T>
void SymmetricTridiagonalSolver<T>::solve_symmetricTridiagonal(T* x, T* scratch) const{

    scratch[0] = sub_diagonal(0) / main_diagonal(0);
    x[0] /= main_diagonal(0);

    for (int i = 1; i < matrix_dimension_-1; i++){
        scratch[i] = sub_diagonal(i) / (main_diagonal(i) - sub_diagonal(i-1) * scratch[i-1]);
        x[i] = (x[i] - sub_diagonal(i-1) * x[i-1]) / (main_diagonal(i) - sub_diagonal(i-1) * scratch[i-1]);
    }

    const int i = matrix_dimension_ - 1;
    x[i] = (x[i] - sub_diagonal(i-1) * x[i-1]) / (main_diagonal(i) - sub_diagonal(i-1) * scratch[i-1]);

    for (int i = matrix_dimension_ - 2; i >= 0; i--){
        x[i] -= scratch[i] * x[i + 1];
    }
}

// ------------------------- //
// Cyclic Tridiagonal Solver //
// ------------------------- //

/* Algorithm based on: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm */
template<typename T>
void SymmetricTridiagonalSolver<T>::solve_symmetricCyclicTridiagonal(T* x, T* u, T* scratch) const{
    const double gamma = -main_diagonal(0);
    const double first_main_diagonal = main_diagonal(0) - gamma;
    const double last_main_diagonal = main_diagonal(matrix_dimension_-1) - cyclic_corner_element() * cyclic_corner_element() / gamma;

    scratch[0] = sub_diagonal(0) / first_main_diagonal;
    x[0] /= first_main_diagonal;
    u[0] = gamma / first_main_diagonal;

    for (int i = 1; i < matrix_dimension_-1; i++){
        const double divisor = 1.0 / (main_diagonal(i) - sub_diagonal(i-1) * scratch[i-1]);
        scratch[i] = sub_diagonal(i) * divisor;
        x[i] = (x[i] - sub_diagonal(i-1) * x[i-1]) * divisor;
        u[i] = (0.0 - sub_diagonal(i-1) * u[i-1]) * divisor;
    }

    const int i = matrix_dimension_ - 1;
    const double divisor = 1.0 / (last_main_diagonal - sub_diagonal(i-1) * scratch[i-1]);
    x[i] = (x[i] - sub_diagonal(i-1) * x[i-1]) * divisor;
    u[i] = (cyclic_corner_element() - sub_diagonal(i-1) * u[i-1]) * divisor;

    for (int i = matrix_dimension_-2; i >= 0; i--){
        x[i] -= scratch[i] * x[i+1];
        u[i] -= scratch[i] * u[i+1];
    }

    const double dot_v_y = x[0] + cyclic_corner_element() / gamma * x[matrix_dimension_-1];
    const double dot_v_q = u[0] + cyclic_corner_element() / gamma * u[matrix_dimension_-1];
    const double factor = dot_v_y / (1.0 + dot_v_q);

    for (int i = 0; i < matrix_dimension_; i++){
        x[i] -= factor * u[i];
    }
}
