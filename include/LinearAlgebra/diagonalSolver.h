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

template<typename T>
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
private:
    int matrix_dimension_;
    std::unique_ptr<T[]> diagonal_values_;
};


// default construction
template<typename T>
DiagonalSolver<T>::DiagonalSolver() :
    matrix_dimension_(0),
    diagonal_values_(nullptr)
{}

// copy construction
template<typename T>
DiagonalSolver<T>::DiagonalSolver(const DiagonalSolver& other) :
    matrix_dimension_(other.matrix_dimension_),
    diagonal_values_(std::make_unique<T[]>(matrix_dimension_))
{
    std::copy(other.diagonal_values_.get(), other.diagonal_values_.get() + matrix_dimension_, diagonal_values_.get());
}

// copy assignment
template<typename T>
DiagonalSolver<T>& DiagonalSolver<T>::operator=(const DiagonalSolver& other){
    if (this == &other) {
        // Self-assignment, no work needed
        return *this;
    }
    // Only allocate new memory if the sizes are different
    if (matrix_dimension_ != other.matrix_dimension_) {
        matrix_dimension_ = other.matrix_dimension_;
        diagonal_values_ = std::make_unique<T[]>(matrix_dimension_);
    }
    std::copy(other.diagonal_values_.get(), other.diagonal_values_.get() + matrix_dimension_, diagonal_values_.get());
    return *this;
}

// move construction
template<typename T>
DiagonalSolver<T>::DiagonalSolver(DiagonalSolver&& other) noexcept :
    matrix_dimension_(other.matrix_dimension_),
    diagonal_values_(std::move(other.diagonal_values_))
{}

// move assignment
template<typename T>
DiagonalSolver<T>& DiagonalSolver<T>::operator=(DiagonalSolver&& other) noexcept{
    matrix_dimension_ = other.matrix_dimension_;
    diagonal_values_ = std::move(other.diagonal_values_);
    return *this;
}

template<typename T>
DiagonalSolver<T>::DiagonalSolver(const int matrix_dimension):
    matrix_dimension_(matrix_dimension),
    diagonal_values_(std::make_unique<T[]>(matrix_dimension_))
{
    assert(matrix_dimension_ >= 1);
}

template<typename T>
int DiagonalSolver<T>::rows() const {
    assert(this->matrix_dimension_ >= 0);
    return this->matrix_dimension_;
}
template<typename T>
int DiagonalSolver<T>::columns() const {
    assert(this->matrix_dimension_ >= 0);
    return this->matrix_dimension_;
}

template<typename T>
const T& DiagonalSolver<T>::diagonal(const int index) const{
    assert(index >= 0); assert(index < this->matrix_dimension_);
    return this->diagonal_values_[index];
}
template<typename T>
T& DiagonalSolver<T>::diagonal(const int index){
    assert(index >= 0); assert(index < this->matrix_dimension_);
    return this->diagonal_values_[index];  
}

// --------------- //
// Diagonal Solver //
// --------------- //

template<typename T>
void DiagonalSolver<T>::solveInPlace(T* sol_rhs) const{
    for (int i = 0; i < matrix_dimension_; i++){
        sol_rhs[i] /= diagonal(i);
    }
}