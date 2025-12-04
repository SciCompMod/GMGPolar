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
class SparseMatrixCOO
{
public:
    using triplet_type = std::tuple<int, int, T>;

    SparseMatrixCOO();
    SparseMatrixCOO(const SparseMatrixCOO& other);
    SparseMatrixCOO(SparseMatrixCOO&& other) noexcept;

    explicit SparseMatrixCOO(int rows, int columns, int nnz);
    explicit SparseMatrixCOO(int rows, int columns, const std::vector<triplet_type>& entries);

    SparseMatrixCOO& operator=(const SparseMatrixCOO& other);
    SparseMatrixCOO& operator=(SparseMatrixCOO&& other) noexcept;

    int rows() const;
    int columns() const;
    int non_zero_size() const;

    const int& row_index(int nz_index) const;
    int& row_index(int nz_index);

    const int& col_index(int nz_index) const;
    int& col_index(int nz_index);

    const T& value(int nz_index) const;
    T& value(int nz_index);

    bool is_symmetric() const;
    void is_symmetric(bool value);

    int* row_indices_data() const;
    int* column_indices_data() const;
    T* values_data() const;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const SparseMatrixCOO<U>& matrix);

    void write_to_file(const std::string& filename) const;

private:
    int rows_;
    int columns_;
    int nnz_;
    AllocatableVector<int> row_indices_;
    AllocatableVector<int> column_indices_;
    AllocatableVector<T> values_;
    bool is_symmetric_ = false;
};

template <typename U>
std::ostream& operator<<(std::ostream& stream, const SparseMatrixCOO<U>& matrix)
{
    stream << "SparseMatrixCOO: " << matrix.rows_ << " x " << matrix.columns_ << "\n";
    stream << "Number of non-zeros (nnz): " << matrix.nnz_ << "\n";
    if (matrix.is_symmetric_) {
        stream << "Matrix is symmetric.\n";
    }
    stream << "Non-zero elements (row, column, value):\n";
    for (int i = 0; i < matrix.nnz_; ++i) {
        stream << "(" << matrix.row_indices_(i) << ", " << matrix.column_indices_(i) << ", " << matrix.values_(i)
               << ")\n";
    }
    return stream;
}

template <typename T>
void SparseMatrixCOO<T>::write_to_file(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");
    }
    file << "SparseMatrixCOO: " << rows_ << " x " << columns_ << "\n";
    file << "Number of non-zeros (nnz): " << nnz_ << "\n";
    if (is_symmetric_) {
        file << "Matrix is symmetric.\n";
    }
    file << "Non-zero elements (row, column, value):\n";
    for (int i = 0; i < nnz_; ++i) {
        file << "(" << row_indices_(i) << ", " << column_indices_(i) << ", " << values_(i) << ")\n";
    }
    file.close();
}

template <typename T>
void sort_entries(std::vector<std::tuple<int, int, T>>& entries)
{
    const auto compare = [](const auto entry1, const auto entry2) {
        const auto local_r1 = std::get<0>(entry1);
        const auto local_r2 = std::get<0>(entry2);
        if (local_r1 < local_r2) {
            return true;
        }
        else if (local_r1 == local_r2) {
            return std::get<1>(entry1) < std::get<1>(entry2);
        }
        return false;
    };
    std::sort(entries.begin(), entries.end(), compare);
}

// default construction
template <typename T>
SparseMatrixCOO<T>::SparseMatrixCOO()
    : rows_(0)
    , columns_(0)
    , nnz_(0)
    , is_symmetric_(false)
{
}

// copy construction
template <typename T>
SparseMatrixCOO<T>::SparseMatrixCOO(const SparseMatrixCOO& other)
    : rows_(other.rows_)
    , columns_(other.columns_)
    , nnz_(other.nnz_)
    , row_indices_("COO row indices", nnz_)
    , column_indices_("COO column indices", nnz_)
    , values_("COO values", nnz_)
    , is_symmetric_(other.is_symmetric_)
{
    Kokkos::deep_copy(row_indices_, other.row_indices_);
    Kokkos::deep_copy(column_indices_, other.column_indices_);
    Kokkos::deep_copy(values_, other.values_);
}

// copy assignment
template <typename T>
SparseMatrixCOO<T>& SparseMatrixCOO<T>::operator=(const SparseMatrixCOO& other)
{
    if (this == &other) {
        // Self-assignment, no work needed
        return *this;
    }
    // Only allocate new memory if the sizes are different
    if (nnz_ != other.nnz_) {
        row_indices_    = Vector<int>("COO row indices", other.nnz_);
        column_indices_ = Vector<int>("COO column indices", other.nnz_);
        values_         = Vector<T>("COO values", other.nnz_);
    }
    // Copy the elements
    rows_         = other.rows_;
    columns_      = other.columns_;
    nnz_          = other.nnz_;
    is_symmetric_ = other.is_symmetric_;
    Kokkos::deep_copy(row_indices_, other.row_indices_);
    Kokkos::deep_copy(column_indices_, other.column_indices_);
    Kokkos::deep_copy(values_, other.values_);
    return *this;
}

// move construction
template <typename T>
SparseMatrixCOO<T>::SparseMatrixCOO(SparseMatrixCOO&& other) noexcept
    : rows_(other.rows_)
    , columns_(other.columns_)
    , nnz_(other.nnz_)
    , row_indices_(std::move(other.row_indices_))
    , column_indices_(std::move(other.column_indices_))
    , values_(std::move(other.values_))
    , is_symmetric_(other.is_symmetric_)
{
    other.nnz_          = 0;
    other.rows_         = 0;
    other.columns_      = 0;
    other.is_symmetric_ = false;
}

// move assignment
template <typename T>
SparseMatrixCOO<T>& SparseMatrixCOO<T>::operator=(SparseMatrixCOO&& other) noexcept
{
    rows_               = other.rows_;
    columns_            = other.columns_;
    nnz_                = other.nnz_;
    row_indices_        = std::move(other.row_indices_);
    column_indices_     = std::move(other.column_indices_);
    values_             = std::move(other.values_);
    is_symmetric_       = other.is_symmetric_;
    other.nnz_          = 0;
    other.rows_         = 0;
    other.columns_      = 0;
    other.is_symmetric_ = false;
    return *this;
}

template <typename T>
SparseMatrixCOO<T>::SparseMatrixCOO(int rows, int columns, int nnz)
    : rows_(rows)
    , columns_(columns)
    , nnz_(nnz)
    , row_indices_("COO row indices", nnz)
    , column_indices_("COO column indices", nnz)
    , values_("COO values", nnz)
    , is_symmetric_(false)
{
    assert(rows >= 0);
    assert(columns >= 0);
    assert(nnz >= 0);
}

template <typename T>
SparseMatrixCOO<T>::SparseMatrixCOO(int rows, int columns, const std::vector<triplet_type>& entries)
    : // entries: row_idx, col_idx, value
    rows_(rows)
    , columns_(columns)
    , nnz_(entries.size())
    , row_indices_("COO row indices", nnz_)
    , column_indices_("COO column indices", nnz_)
    , values_("COO values", nnz_)
    , is_symmetric_(false)
{
    assert(rows_ >= 0);
    assert(columns_ >= 0);
    assert(nnz_ >= 0);
#pragma omp parallel for
    for (int i = 0; i < nnz_; i++) {
        assert(0 <= std::get<0>(entries[i]) && std::get<0>(entries[i]) < rows_);
        assert(0 <= std::get<1>(entries[i]) && std::get<1>(entries[i]) < columns_);
        row_indices_(i)    = std::get<0>(entries[i]);
        column_indices_(i) = std::get<1>(entries[i]);
        values_(i)         = std::get<2>(entries[i]);
    }
}

template <typename T>
int SparseMatrixCOO<T>::rows() const
{
    assert(this->rows_ >= 0);
    return this->rows_;
}
template <typename T>
int SparseMatrixCOO<T>::columns() const
{
    assert(this->columns_ >= 0);
    return this->columns_;
}
template <typename T>
int SparseMatrixCOO<T>::non_zero_size() const
{
    assert(this->nnz_ >= 0);
    assert(static_cast<size_t>(this->nnz_) <= static_cast<size_t>(this->rows_) * static_cast<size_t>(this->columns_));
    return this->nnz_;
}

template <typename T>
int& SparseMatrixCOO<T>::row_index(int nz_index)
{
    assert(nz_index >= 0);
    assert(nz_index < this->nnz_);
    return row_indices_(nz_index);
}
template <typename T>
const int& SparseMatrixCOO<T>::row_index(int nz_index) const
{
    assert(nz_index >= 0);
    assert(nz_index < this->nnz_);
    return row_indices_(nz_index);
}

template <typename T>
int& SparseMatrixCOO<T>::col_index(int nz_index)
{
    assert(nz_index >= 0);
    assert(nz_index < nnz_);
    return column_indices_(nz_index);
}
template <typename T>
const int& SparseMatrixCOO<T>::col_index(int nz_index) const
{
    assert(nz_index >= 0);
    assert(nz_index < nnz_);
    return column_indices_(nz_index);
}

template <typename T>
T& SparseMatrixCOO<T>::value(int nz_index)
{
    assert(nz_index >= 0);
    assert(nz_index < nnz_);
    return values_(nz_index);
}
template <typename T>
const T& SparseMatrixCOO<T>::value(int nz_index) const
{
    assert(nz_index >= 0);
    assert(nz_index < nnz_);
    return values_(nz_index);
}

template <typename T>
bool SparseMatrixCOO<T>::is_symmetric() const
{
    return is_symmetric_;
}
template <typename T>
void SparseMatrixCOO<T>::is_symmetric(bool value)
{
    is_symmetric_ = value;
}

template <typename T>
int* SparseMatrixCOO<T>::row_indices_data() const
{
    return row_indices_.data();
}

template <typename T>
int* SparseMatrixCOO<T>::column_indices_data() const
{
    return column_indices_.data();
}

template <typename T>
T* SparseMatrixCOO<T>::values_data() const
{
    return values_.data();
}
