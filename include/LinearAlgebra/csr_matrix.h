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

// The CSR matrix format is used if a custom LU decomposition solver is choosen.
// MUMPS relies on the COO format.

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
    explicit SparseMatrixCSR(int rows, int columns, const std::vector<T>& values,
                             const std::vector<int>& column_indices, const std::vector<int>& row_start_indices);

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
    AllocatableVector<T> values_;
    AllocatableVector<int> column_indices_;
    AllocatableVector<int> row_start_indices_;

    bool is_sorted_entries(const std::vector<std::tuple<int, int, T>>& entries)
    {
        for (size_t i = 1; i < entries.size(); ++i) {
            const auto& prev = entries[i - 1];
            const auto& curr = entries[i];
            if (std::get<0>(prev) > std::get<0>(curr))
                return false;
        }
        return true;
    }
};

template <typename U>
std::ostream& operator<<(std::ostream& stream, const SparseMatrixCSR<U>& matrix)
{
    stream << "SparseMatrixCSR: " << matrix.rows_ << " x " << matrix.columns_ << "\n";
    stream << "Number of non-zeros (nnz): " << matrix.nnz_ << "\n";
    stream << "Non-zero elements (row, column, value):\n";
    for (int row = 0; row < matrix.rows_; ++row) {
        for (int nnz = matrix.row_start_indices_(row); nnz < matrix.row_start_indices_(row + 1); ++nnz) {
            stream << "(" << row << ", " << matrix.column_indices_(nnz) << ", " << matrix.values_(nnz) << ")\n";
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
{
}

// copy construction
template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(const SparseMatrixCSR& other)
    : rows_(other.rows_)
    , columns_(other.columns_)
    , nnz_(other.nnz_)
    , values_("CSR values", nnz_)
    , column_indices_("CSR column indices", nnz_)
    , row_start_indices_("CSR row start indices", rows_ + 1)
{
    Kokkos::deep_copy(values_, other.values_);
    Kokkos::deep_copy(column_indices_, other.column_indices_);
    Kokkos::deep_copy(row_start_indices_, other.row_start_indices_);
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
        values_            = Vector<T>("CSR values", other.nnz_);
        column_indices_    = Vector<int>("CSR column indices", other.nnz_);
        row_start_indices_ = Vector<int>("CSR row start indices", other.rows_ + 1);
    }
    // Copy the elements
    rows_    = other.rows_;
    columns_ = other.columns_;
    nnz_     = other.nnz_;
    Kokkos::deep_copy(values_, other.values_);
    Kokkos::deep_copy(column_indices_, other.column_indices_);
    Kokkos::deep_copy(row_start_indices_, other.row_start_indices_);
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
    rows_              = other.rows_;
    columns_           = other.columns_;
    nnz_               = other.nnz_;
    values_            = std::move(other.values_);
    column_indices_    = std::move(other.column_indices_);
    row_start_indices_ = std::move(other.row_start_indices_);
    other.nnz_         = 0;
    other.rows_        = 0;
    other.columns_     = 0;
    return *this;
}

template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(int rows, int columns, std::function<int(int)> nz_per_row)
    : rows_(rows)
    , columns_(columns)
    , row_start_indices_("CSR row start indices", rows_ + 1)
{
    assert(rows >= 0);
    assert(columns >= 0);
    nnz_ = 0;
    for (int i = 0; i < rows; i++) {
        row_start_indices_(i) = nnz_;
        nnz_ += nz_per_row(i);
    }
    row_start_indices_(rows) = nnz_;
    values_                  = Vector<T>("CSR values", nnz_);
    column_indices_          = Vector<int>("CSR column indices", nnz_);
}

template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(int rows, int columns, const std::vector<triplet_type>& entries)
    : // entries: row_idx, col_idx, value
    rows_(rows)
    , columns_(columns)
    , nnz_(entries.size())
    , values_("CSR values", nnz_)
    , column_indices_("CSR column indices", nnz_)
    , row_start_indices_("CSR row start indices", rows_ + 1)
{
    assert(rows >= 0);
    assert(columns >= 0);
    assert(is_sorted_entries(entries) && "Entries must be sorted by row!");
    // fill values and column indexes
    for (int i = 0; i < nnz_; i++) {
        assert(0 <= std::get<0>(entries[i]) && std::get<0>(entries[i]) < rows);
        values_(i)         = std::get<2>(entries[i]);
        column_indices_(i) = std::get<1>(entries[i]);
        assert(0 <= column_indices_(i) && column_indices_(i) < columns);
    }
    //fill row indexes
    int count             = 0;
    row_start_indices_(0) = 0;
    for (int r = 0; r < rows; r++) {
        while (count < nnz_ && std::get<0>(entries[count]) == r)
            count++;
        row_start_indices_(r + 1) = count;
    }
    assert(row_start_indices_(rows) == nnz_);
}

template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(int rows, int columns, const std::vector<T>& values,
                                    const std::vector<int>& column_indices, const std::vector<int>& row_start_indices)
    : rows_(rows)
    , columns_(columns)
    , nnz_(values.size())
    , values_("CSR values", nnz_)
    , column_indices_("CSR column indices", nnz_)
    , row_start_indices_("CSR row start indices", rows_ + 1)
{
    assert(rows >= 0);
    assert(columns >= 0);
    assert(row_start_indices.size() == static_cast<size_t>(rows + 1));
    assert(values.size() == column_indices.size());

    // Copy data to internal storage
    std::copy(values.begin(), values.end(), values_.data());
    std::copy(column_indices.begin(), column_indices.end(), column_indices_.data());
    std::copy(row_start_indices.begin(), row_start_indices.end(), row_start_indices_.data());
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
    return row_start_indices_(row + 1) - row_start_indices_(row);
}

template <typename T>
const int& SparseMatrixCSR<T>::row_nz_index(int row, int nz_index) const
{
    assert(row >= 0 && row < rows_);
    assert(nz_index >= 0 && nz_index < row_nz_size(row));
    return column_indices_(row_start_indices_(row) + nz_index);
}

template <typename T>
int& SparseMatrixCSR<T>::row_nz_index(int row, int nz_index)
{
    assert(row >= 0 && row < rows_);
    assert(nz_index >= 0 && nz_index < row_nz_size(row));
    return column_indices_(row_start_indices_(row) + nz_index);
}

template <typename T>
const T& SparseMatrixCSR<T>::row_nz_entry(int row, int nz_index) const
{
    assert(row >= 0 && row < rows_);
    assert(nz_index >= 0 && nz_index < row_nz_size(row));
    return values_(row_start_indices_(row) + nz_index);
}

template <typename T>
T& SparseMatrixCSR<T>::row_nz_entry(int row, int nz_index)
{
    assert(row >= 0 && row < rows_);
    assert(nz_index >= 0 && nz_index < row_nz_size(row));
    return values_(row_start_indices_(row) + nz_index);
}

template <typename T>
T* SparseMatrixCSR<T>::values_data() const
{
    return values_.data();
}

template <typename T>
int* SparseMatrixCSR<T>::column_indices_data() const
{
    return column_indices_.data();
}

template <typename T>
int* SparseMatrixCSR<T>::row_start_indices_data() const
{
    return row_start_indices_.data();
}
