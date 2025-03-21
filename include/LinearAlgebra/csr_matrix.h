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
#include <vector>
#include <unordered_map>
#include <cassert>
#include <iostream>
#include <cmath>

#include <fstream>
#include <iostream>

#include "vector.h"

/* The CSR matrix format is currently unused, as we use MUMPS which relies on the COO format. */
/* Here we provide a custom LU decomposition solver, which could be replaced by different library implementation, 
if we would decide to move away from mumps. */

template <typename T>
class SparseMatrixCSR
{
public:
    using triplet_type = std::tuple<int, int, T>;

    SparseMatrixCSR();
    SparseMatrixCSR(const SparseMatrixCSR& other);
    SparseMatrixCSR(SparseMatrixCSR&& other) noexcept;

    explicit SparseMatrixCSR(int rows, int columns, std::function<int(int)> nz_per_row);
    explicit SparseMatrixCSR(int rows, int columns, const std::vector<triplet_type>& entries);
    explicit SparseMatrixCSR(int rows, int columns, 
        const std::vector<T>& values, const std::vector<int>& column_indices, const std::vector<int>& row_start_indices
    );

    SparseMatrixCSR& operator=(const SparseMatrixCSR& other);
    SparseMatrixCSR& operator=(SparseMatrixCSR&& other) noexcept;

    int rows() const;
    int columns() const;
    int non_zero_size() const;

    int row_nz_size(int row) const;

    const int& row_nz_index(int row, int nz_index) const;
    int& row_nz_index(int row, int nz_index);

    const T& row_nz_entry(int row, int nz_index) const;
    T& row_nz_entry(int row, int nz_index);

    T* values_data() const;
    int* column_indices_data() const;
    int* row_start_indices_data() const;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const SparseMatrixCSR<U>& matrix);

private:
    int rows_;
    int columns_;
    int nnz_;
    std::unique_ptr<T[]> values_;
    std::unique_ptr<int[]> column_indices_;
    std::unique_ptr<int[]> row_start_indices_;
};

/* LU decomposition Solver (slow) */
template <typename T>
class SparseLUSolver {
public:
    SparseLUSolver();
    SparseLUSolver(const SparseMatrixCSR<T>& A);
    void solveInPlace(Vector<T>& b) const;
private:
    std::vector<T> L_values, U_values;
    std::vector<int> L_col_idx, U_col_idx;
    std::vector<int> L_row_ptr, U_row_ptr;
    bool factorized_ = false;
    void factorize(const SparseMatrixCSR<T>& A);
};







template <typename U>
std::ostream& operator<<(std::ostream& stream, const SparseMatrixCSR<U>& matrix)
{
    stream << "SparseMatrixCSR: " << matrix.rows_ << " x " << matrix.columns_ << "\n";
    stream << "Number of non-zeros (nnz): " << matrix.nnz_ << "\n";
    stream << "Non-zero elements (row, column, value):\n";
    for (int row = 0; row < matrix.rows_; ++row) {
        for (int nnz = matrix.row_start_indices_[row]; nnz < matrix.row_start_indices_[row + 1]; ++nnz) {
            stream << "(" << row << ", " << matrix.column_indices_[nnz] << ", " << matrix.values_[nnz] << ")\n";
        }
    }
    return stream;
}

// default construction
template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR()
    : rows_(0)
    , columns_(0)
    , nnz_(0)
    , values_(nullptr)
    , column_indices_(nullptr)
    , row_start_indices_(nullptr)
{
}

// copy construction
template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(const SparseMatrixCSR& other)
    : rows_(other.rows_)
    , columns_(other.columns_)
    , nnz_(other.nnz_)
    , values_(std::make_unique<T[]>(nnz_))
    , column_indices_(std::make_unique<int[]>(nnz_))
    , row_start_indices_(std::make_unique<int[]>(rows_ + 1))
{
    std::copy(other.values_.get(), other.values_.get() + nnz_, values_.get());
    std::copy(other.column_indices_.get(), other.column_indices_.get() + nnz_, column_indices_.get());
    std::copy(other.row_start_indices_.get(), other.row_start_indices_.get() + rows_ + 1, row_start_indices_.get());
}

// copy assignment
template <typename T>
SparseMatrixCSR<T>& SparseMatrixCSR<T>::operator=(const SparseMatrixCSR& other)
{
    if (this == &other) {
        // Self-assignment, no work needed
        return *this;
    }
    // Only allocate new memory if the sizes are different
    if (nnz_ != other.nnz_ || rows_ != other.rows_) {
        values_            = std::make_unique<T[]>(other.nnz_);
        column_indices_    = std::make_unique<int[]>(other.nnz_);
        row_start_indices_ = std::make_unique<int[]>(other.rows_ + 1);
    }
    // Copy the elements
    rows_    = other.rows_;
    columns_ = other.columns_;
    nnz_     = other.nnz_;
    std::copy(other.values_.get(), other.values_.get() + nnz_, values_.get());
    std::copy(other.column_indices_.get(), other.column_indices_.get() + nnz_, column_indices_.get());
    std::copy(other.row_start_indices_.get(), other.row_start_indices_.get() + rows_ + 1, row_start_indices_.get());
    return *this;
}

// move construction
template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(SparseMatrixCSR&& other) noexcept
    : rows_(other.rows_)
    , columns_(other.columns_)
    , nnz_(other.nnz_)
    , values_(std::move(other.values_))
    , column_indices_(std::move(other.column_indices_))
    , row_start_indices_(std::move(other.row_start_indices_))
{
    other.nnz_     = 0;
    other.rows_    = 0;
    other.columns_ = 0;
}

// move assignment
template <typename T>
SparseMatrixCSR<T>& SparseMatrixCSR<T>::operator=(SparseMatrixCSR&& other) noexcept
{
    rows_               = other.rows_;
    columns_            = other.columns_;
    nnz_                = other.nnz_;
    values_             = std::move(other.values_);
    column_indices_     = std::move(other.column_indices_);
    row_start_indices_  = std::move(other.row_start_indices_);
    other.nnz_          = 0;
    other.rows_         = 0;
    other.columns_      = 0;
    return *this;
}

template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(int rows, int columns, std::function<int(int)> nz_per_row)
    : rows_(rows)
    , columns_(columns)
    , row_start_indices_(std::make_unique<int[]>(rows_ + 1))
{
    assert(rows >= 0);
    assert(columns >= 0);
    nnz_ = 0;
    for (int i = 0; i < rows; i++){
        row_start_indices_[i] = nnz_;
        nnz_ += nz_per_row(i);
    }
    row_start_indices_[rows] = nnz_;
    values_ = std::make_unique<T[]>(nnz_);
    column_indices_ = std::make_unique<int[]>(nnz_);
}

template<typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(int rows, int columns, const std::vector<triplet_type>& entries):
    // entries: row_idx, col_idx, value
    rows_(rows),
    columns_(columns),
    nnz_(entries.size()),
    values_(std::make_unique<T[]>(nnz_)),
    column_indices_(std::make_unique<int[]>(nnz_)),
    row_start_indices_(std::make_unique<int[]>(rows_ + 1))
{
    assert(rows >= 0);
    assert(columns >= 0);
    // fill values and column indexes
    for (int i = 0; i < nnz_; i++){
        assert(0 <= std::get<0>(entries[i]) && std::get<0>(entries[i]) < rows);
        values_[i] = std::get<2>(entries[i]);
        column_indices_[i] = std::get<1>(entries[i]);
        assert(0 <= column_indices_[i] && column_indices_[i] < columns);
    }
    //fill row indexes
    int count = 0;
    row_start_indices_[0] = 0;
    for (int r = 0; r < rows; r++){
        while (count < nnz_ && std::get<0>(entries[count]) == r) count++;
        row_start_indices_[r + 1] = count;
    }
    assert(row_start_indices_[rows] == nnz_);
}

template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(
    int rows, int columns, 
    const std::vector<T>& values, 
    const std::vector<int>& column_indices, 
    const std::vector<int>& row_start_indices
) 
    : rows_(rows), columns_(columns), nnz_(values.size()), 
      values_(std::make_unique<T[]>(nnz_)), 
      column_indices_(std::make_unique<int[]>(nnz_)), 
      row_start_indices_(std::make_unique<int[]>(rows_ + 1)) 
{
    assert(rows >= 0);
    assert(columns >= 0);
    assert(row_start_indices.size() == static_cast<size_t>(rows + 1));
    assert(values.size() == column_indices.size());

    // Copy data to internal storage
    std::copy(values.begin(), values.end(), values_.get());
    std::copy(column_indices.begin(), column_indices.end(), column_indices_.get());
    std::copy(row_start_indices.begin(), row_start_indices.end(), row_start_indices_.get());
}


template <typename T>
int SparseMatrixCSR<T>::rows() const
{
    assert(this->rows_ >= 0);
    return this->rows_;
}
template <typename T>
int SparseMatrixCSR<T>::columns() const
{
    assert(this->columns_ >= 0);
    return this->columns_;
}
template <typename T>
int SparseMatrixCSR<T>::non_zero_size() const
{
    assert(this->nnz_ >= 0);
    assert(static_cast<size_t>(this->nnz_) <= static_cast<size_t>(this->rows_) * static_cast<size_t>(this->columns_));
    return this->nnz_;
}

template <typename T>
int SparseMatrixCSR<T>::row_nz_size(int row) const
{
    assert(row >= 0 && row < rows_);
    return row_start_indices_[row + 1] - row_start_indices_[row];
}

template <typename T>
const int& SparseMatrixCSR<T>::row_nz_index(int row, int nz_index) const {
    assert(row >= 0 && row < rows_);
    assert(nz_index >= 0 && nz_index < row_nz_size(row));
    return column_indices_[row_start_indices_[row] + nz_index];
}

template <typename T>
int& SparseMatrixCSR<T>::row_nz_index(int row, int nz_index) {
    assert(row >= 0 && row < rows_);
    assert(nz_index >= 0 && nz_index < row_nz_size(row));
    return column_indices_[row_start_indices_[row] + nz_index];
}

template <typename T>
const T& SparseMatrixCSR<T>::row_nz_entry(int row, int nz_index) const {
    assert(row >= 0 && row < rows_);
    assert(nz_index >= 0 && nz_index < row_nz_size(row));
    return values_[row_start_indices_[row] + nz_index];
}

template <typename T>
T& SparseMatrixCSR<T>::row_nz_entry(int row, int nz_index) {
    assert(row >= 0 && row < rows_);
    assert(nz_index >= 0 && nz_index < row_nz_size(row));
    return values_[row_start_indices_[row] + nz_index];
}

template <typename T>
T* SparseMatrixCSR<T>::values_data() const {
    return values_.get();
}

template <typename T>
int* SparseMatrixCSR<T>::column_indices_data() const {
    return column_indices_.get();
}

template <typename T>
int* SparseMatrixCSR<T>::row_start_indices_data() const {
    return row_start_indices_.get();
}



template <typename T>
SparseLUSolver<T>::SparseLUSolver() 
    : factorized_(false) {
}

template <typename T>
SparseLUSolver<T>::SparseLUSolver(const SparseMatrixCSR<T>& A) {
    assert(A.rows() == A.columns());
    factorize(A);
}


/* This is slow */
template <typename T>
void SparseLUSolver<T>::factorize(const SparseMatrixCSR<T>& A) {
    const int n = A.rows();
    L_row_ptr.resize(n + 1, 0);
    U_row_ptr.resize(n + 1, 0);

    // Temporary structures to store computed rows
    std::vector<std::unordered_map<int, T>> L_map(n);
    std::vector<std::unordered_map<int, T>> U_map(n);
    
    for (int i = 0; i < n; i++) {
        std::unordered_map<int, T> row_values;
        
        // Load nonzero elements of row i
        for (int idx = 0; idx < A.row_nz_size(i); idx++) {
            int j = A.row_nz_index(i, idx);
            row_values[j] = A.row_nz_entry(i, idx);
        }

        for (int j = 0; j < i; j++) {
            if (row_values.find(j) == row_values.end()) continue;

            row_values[j] /= U_map[j][j];
            L_map[i][j] = row_values[j];

            for (const auto& [k, U_val] : U_map[j]) {
                if (k > j) {
                    row_values[k] -= row_values[j] * U_val;
                }
            }
        }

        for (const auto& [j, val] : row_values) {
            if (j < i)
                L_map[i][j] = val;
            else
                U_map[i][j] = val;
        }
    }

    // Convert L_map and U_map into CSR format
    for (int i = 0; i < n; i++) {
        for (const auto& [col, val] : L_map[i]) {
            L_values.push_back(val);
            L_col_idx.push_back(col);
        }
        L_row_ptr[i + 1] = L_values.size();

        for (const auto& [col, val] : U_map[i]) {
            U_values.push_back(val);
            U_col_idx.push_back(col);
        }
        U_row_ptr[i + 1] = U_values.size();
    }

    factorized_ = true;
}

template <typename T>
void SparseLUSolver<T>::solveInPlace(Vector<T>& b) const {
    assert(factorized_);
    const int n = L_row_ptr.size() - 1;
    assert(b.size() == n);

    // Forward substitution (L * b = b) -> b now holds y
    for (int i = 0; i < n; i++) {
        for (int idx = L_row_ptr[i]; idx < L_row_ptr[i + 1]; idx++) {
            b[i] -= L_values[idx] * b[L_col_idx[idx]];
        }
    }
    
    // Backward substitution (U * b = y) -> b now holds x
    for (int i = n - 1; i >= 0; i--) {
        T diag = 0;
        for (int idx = U_row_ptr[i]; idx < U_row_ptr[i + 1]; idx++) {
            int col = U_col_idx[idx];
            if (col == i) {
                diag = U_values[idx];  // Store diagonal value
            } else {
                b[i] -= U_values[idx] * b[col];
            }
        }
        if (std::abs(diag) < 1e-12) {
            std::cerr << "Zero diagonal encountered in U at row " << i << "!\n";
            std::exit(EXIT_FAILURE);
        }
        b[i] /= diag;
    }
}
