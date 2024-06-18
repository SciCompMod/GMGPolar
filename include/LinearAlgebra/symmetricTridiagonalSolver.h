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
    void SolveInPlace(T* sol_rhs, T* temp1, T* temp2 = nullptr);
private:
    int matrix_dimension_;
    std::unique_ptr<T[]> main_diagonal_values_;
    std::unique_ptr<T[]> sub_diagonal_values_;
    T cyclic_corner_element_ = 0.0;
    bool is_cyclic_ = true;
    bool solver_setup_ = false;

    // Solve methods
    void solve_symmetricTridiagonal(T* x, T* scratch);
    void solve_symmetricCyclicTridiagonal(T* x, T* chmod, T* u);
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
    matrix_dimension_ = other.matrix_dimension_;
    main_diagonal_values_ = std::make_unique<T[]>(matrix_dimension_);
    sub_diagonal_values_ = std::make_unique<T[]>(matrix_dimension_-1);
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
    return this->cyclic_corner_element_;
}
template<typename T>
T& SymmetricTridiagonalSolver<T>::cyclic_corner_element(){
    return this->cyclic_corner_element_;
}

template<typename T>
void SymmetricTridiagonalSolver<T>::SolveInPlace(T* sol_rhs, T* temp1, T* temp2){
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
void SymmetricTridiagonalSolver<T>::solve_symmetricTridiagonal(T* x, T* scratch){
    assert(!equals(main_diagonal(0), 0.0));

    scratch[0] = sub_diagonal(0) / main_diagonal(0);
    x[0] /= main_diagonal(0);

    for (int i = 1; i < matrix_dimension_; i++){
        if(i < matrix_dimension_ - 1){
            assert(!equals(main_diagonal(i) - sub_diagonal(i-1) * scratch[i-1], 0.0));
            scratch[i] = sub_diagonal(i) / (main_diagonal(i) - sub_diagonal(i-1) * scratch[i-1]);
        }
        assert(!equals(main_diagonal(i) - sub_diagonal(i-1) * scratch[i-1], 0.0));
        x[i] = (x[i] - sub_diagonal(i-1) * x[i-1]) / (main_diagonal(i) - sub_diagonal(i-1) * scratch[i-1]);
    }

    for (int i = matrix_dimension_ - 2; i >= 0; i--){
        x[i] -= scratch[i] * x[i + 1];
    }
}

// ------------------------- //
// Cyclic Tridiagonal Solver //
// ------------------------- //

/* Algorithm based on: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm */
template<typename T>
void SymmetricTridiagonalSolver<T>::solve_symmetricCyclicTridiagonal(T* x, T* cmod, T* u){
    /* ------------------------------------- */
    /* Optimized version found on wikipedia. */
    /* ------------------------------------- */
    
    /* arbitrary, but chosen such that division by zero is avoided */
    const double gamma = -main_diagonal(0);

    cmod[0] = cyclic_corner_element() / (main_diagonal(0) - gamma);
    u[0] = gamma / (main_diagonal(0) - gamma);
    x[0] /= (main_diagonal(0) - gamma);

    for (int i = 1; i + 1 < matrix_dimension_; i++){
        const double m = 1.0 / (main_diagonal(i) - sub_diagonal(i-1) * cmod[i-1]);
        cmod[i] = sub_diagonal(i) * m;
        u[i] = (0.0 - sub_diagonal(i-1) * u[i-1]) * m;
        x[i] = (x[i] - sub_diagonal(i-1) * x[i-1]) * m;
    }

    const double m = 1.0 / (main_diagonal(matrix_dimension_-1) - cyclic_corner_element() * cyclic_corner_element() / gamma - cyclic_corner_element() * cmod[matrix_dimension_-2]);
    u[matrix_dimension_-1] = (cyclic_corner_element() - sub_diagonal(matrix_dimension_-2) * u[matrix_dimension_-2]) * m;
    x[matrix_dimension_-1] = (x[matrix_dimension_-1] - sub_diagonal(matrix_dimension_-2) * x[matrix_dimension_-2]) * m;

    for (int i = matrix_dimension_-2; i >= 0; i--){
        u[i] -= cmod[i] * u[i+1];
        x[i] -= cmod[i] * x[i+1];
    }

    const double fact = (x[0] + x[matrix_dimension_-1] * cyclic_corner_element() / gamma) / (1.0 + u[0] + u[matrix_dimension_-1] * cyclic_corner_element() / gamma);
    
    for (int i = 0; i < matrix_dimension_; i++){
        x[i] -= fact * u[i];
    }

    /* ------------------------------------------------------------------- */
    /* This is the straigt forward unoptimized version found on wikipedia. */
    /* ------------------------------------------------------------------- */

    /*
    const double gamma = -main_diagonal(0);

    if(!solver_setup_){
        main_diagonal(0) -= gamma;
        main_diagonal(matrix_dimension_-1) -= cyclic_corner_element() * cyclic_corner_element() / gamma;
        solver_setup_ = true;
    }

    solve_symmetricTridiagonal(x, cmod);

    u[0] = gamma;
    for (int i = 1; i < matrix_dimension_-1; i++) u[i] = 0.0;
    u[matrix_dimension_-1] = cyclic_corner_element();

    solve_symmetricTridiagonal(u, cmod);

    const double dot_v_y = x[0] + cyclic_corner_element()/gamma * x[matrix_dimension_-1];
    const double dot_v_q = u[0] + cyclic_corner_element()/gamma * u[matrix_dimension_-1];

    for (int i = 0; i < matrix_dimension_; i++){
        x[i] -= u[i] * dot_v_y / dot_v_q;
    }
    */
}