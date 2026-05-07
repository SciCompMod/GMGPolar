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
#include <unordered_map>
#include <cassert>
#include <iostream>
#include <cmath>

#include <fstream>
#include <iostream>

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

// The CSR matrix format is used if a custom LU decomposition solver is choosen.
// MUMPS relies on the COO format.

namespace gmgpolar
{

template <typename T>
class SparseMatrixCSR
{
public:
    using non_const_element_type = std::remove_const_t<T>;
    using int_element_type       = std::conditional_t<std::is_const_v<T>, const int, int>;

    using const_type     = SparseMatrixCSR<const T>;
    using non_const_type = SparseMatrixCSR<non_const_element_type>;

    using triplet_type = std::tuple<int, int, non_const_element_type>;

    SparseMatrixCSR();
    KOKKOS_DEFAULTED_FUNCTION SparseMatrixCSR(const SparseMatrixCSR& other) = default;
    KOKKOS_DEFAULTED_FUNCTION SparseMatrixCSR(SparseMatrixCSR&& other) noexcept;
    template <typename T2,
              std::enable_if_t<std::is_same_v<non_const_element_type, T2> && std::is_const_v<T>, bool> = true>
    KOKKOS_FUNCTION SparseMatrixCSR(const SparseMatrixCSR<T2>& other);

    explicit SparseMatrixCSR(int rows, int columns, std::function<int(int)> nz_per_row);
    explicit SparseMatrixCSR(int rows, int columns, const std::vector<triplet_type>& entries);
    explicit SparseMatrixCSR(int rows, int columns, const std::vector<non_const_element_type>& values,
                             const std::vector<int>& column_indices, const std::vector<int>& row_start_indices);

    SparseMatrixCSR& operator=(const SparseMatrixCSR& other)                               = delete;
    KOKKOS_DEFAULTED_FUNCTION SparseMatrixCSR& operator=(SparseMatrixCSR&& other) noexcept = default;

    template <typename T2,
              std::enable_if_t<std::is_same_v<std::remove_const_t<T>, T2> && std::is_const_v<T>, bool> = true>
    KOKKOS_FUNCTION SparseMatrixCSR& operator=(const SparseMatrixCSR<T2>& other);

    non_const_type copy() const;
    const_type to_const() const;

    KOKKOS_FUNCTION int rows() const;
    KOKKOS_FUNCTION int columns() const;
    KOKKOS_FUNCTION int non_zero_size() const;

    KOKKOS_FUNCTION int row_nz_size(int row) const;

    KOKKOS_FUNCTION int_element_type& row_nz_index(int row, int nz_index) const;

    KOKKOS_FUNCTION T& row_nz_entry(int row, int nz_index) const;

    KOKKOS_FUNCTION T* values_data() const;
    KOKKOS_FUNCTION int_element_type* column_indices_data() const;
    KOKKOS_FUNCTION int_element_type* row_start_indices_data() const;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const SparseMatrixCSR<U>& matrix);

    template <typename T2>
    friend class SparseMatrixCSR;

private:
    int rows_;
    int columns_;
    int nnz_;
    AllocatableVector<T> values_;
    AllocatableVector<int_element_type> column_indices_;
    AllocatableVector<int_element_type> row_start_indices_;

    bool is_sorted_entries(const std::vector<std::tuple<int, int, non_const_element_type>>& entries)
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
SparseMatrixCSR<T>::non_const_type SparseMatrixCSR<T>::copy() const
{
    non_const_type new_copy;
    new_copy.rows_              = rows_;
    new_copy.columns_           = columns_;
    new_copy.nnz_               = nnz_;
    new_copy.values_            = AllocatableVector<non_const_element_type>("CSR values", nnz_);
    new_copy.column_indices_    = AllocatableVector<int>("CSR column indices", nnz_);
    new_copy.row_start_indices_ = AllocatableVector<int>("CSR row start indices", rows_ + 1);
    Kokkos::deep_copy(new_copy.values_, values_);
    Kokkos::deep_copy(new_copy.column_indices_, column_indices_);
    Kokkos::deep_copy(new_copy.row_start_indices_, row_start_indices_);

    return new_copy;
}

// move construction
template <typename T>
KOKKOS_FUNCTION SparseMatrixCSR<T>::SparseMatrixCSR(SparseMatrixCSR&& other) noexcept
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

// const cast construction
template <typename T>
template <typename T2, std::enable_if_t<std::is_same_v<std::remove_const_t<T>, T2> && std::is_const_v<T>, bool>>
KOKKOS_FUNCTION SparseMatrixCSR<T>::SparseMatrixCSR(const SparseMatrixCSR<T2>& other)
    : rows_(other.rows_)
    , columns_(other.columns_)
    , nnz_(other.nnz_)
    , values_(other.values_)
    , column_indices_(other.column_indices_)
    , row_start_indices_(other.row_start_indices_)
{
}

template <typename T>
template <typename T2, std::enable_if_t<std::is_same_v<std::remove_const_t<T>, T2> && std::is_const_v<T>, bool>>
KOKKOS_FUNCTION SparseMatrixCSR<T>& SparseMatrixCSR<T>::operator=(const SparseMatrixCSR<T2>& other)
{
    rows_              = other.rows_;
    columns_           = other.columns_;
    nnz_               = other.nnz_;
    values_            = other.values_;
    column_indices_    = other.column_indices_;
    row_start_indices_ = other.row_start_indices_;
    return *this;
}

template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(int rows, int columns, std::function<int(int)> nz_per_row)
    : rows_(rows)
    , columns_(columns)
{
    assert(rows >= 0);
    assert(columns >= 0);

    Vector<int> row_start_indices("CSR row start indices", rows_ + 1);

    nnz_ = 0;
    for (int i = 0; i < rows; i++) {
        row_start_indices(i) = nnz_;
        nnz_ += nz_per_row(i);
    }
    row_start_indices(rows) = nnz_;

    Vector<non_const_element_type> values("CSR values", nnz_);
    Vector<int> column_indices("CSR column indices", nnz_);

    assign(values, T(0));
    assign(column_indices, 0);

    values_            = values;
    column_indices_    = column_indices;
    row_start_indices_ = row_start_indices;
}

template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(int rows, int columns, const std::vector<triplet_type>& entries)
    : // entries: row_idx, col_idx, value
    rows_(rows)
    , columns_(columns)
    , nnz_(entries.size())
{
    assert(rows >= 0);
    assert(columns >= 0);
    assert(is_sorted_entries(entries) && "Entries must be sorted by row!");

    Vector<non_const_element_type> values("CSR values", nnz_);
    Vector<int> column_indices("CSR column indices", nnz_);
    Vector<int> row_start_indices("CSR row start indices", rows_ + 1);

    // fill values and column indexes
    for (int i = 0; i < nnz_; i++) {
        assert(0 <= std::get<0>(entries[i]) && std::get<0>(entries[i]) < rows);
        values(i)         = std::get<2>(entries[i]);
        column_indices(i) = std::get<1>(entries[i]);
        assert(0 <= column_indices(i) && column_indices(i) < columns);
    }
    //fill row indexes
    int count            = 0;
    row_start_indices(0) = 0;
    for (int r = 0; r < rows; r++) {
        while (count < nnz_ && std::get<0>(entries[count]) == r)
            count++;
        row_start_indices(r + 1) = count;
    }
    assert(row_start_indices(rows) == nnz_);

    values_            = values;
    column_indices_    = column_indices;
    row_start_indices_ = row_start_indices;
}

template <typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(int rows, int columns, const std::vector<non_const_element_type>& values,
                                    const std::vector<int>& column_indices, const std::vector<int>& row_start_indices)
    : rows_(rows)
    , columns_(columns)
    , nnz_(values.size())
{
    assert(rows >= 0);
    assert(columns >= 0);
    assert(row_start_indices.size() == static_cast<size_t>(rows + 1));
    assert(values.size() == column_indices.size());

    Vector<non_const_element_type> csr_values("CSR values", nnz_);
    Vector<int> csr_column_indices("CSR column indices", nnz_);
    Vector<int> csr_row_start_indices("CSR row start indices", rows_ + 1);

    // Copy data to internal storage
    std::copy(values.begin(), values.end(), csr_values.data());
    std::copy(column_indices.begin(), column_indices.end(), csr_column_indices.data());
    std::copy(row_start_indices.begin(), row_start_indices.end(), csr_row_start_indices.data());

    values_            = csr_values;
    column_indices_    = csr_column_indices;
    row_start_indices_ = csr_row_start_indices;
}

template <typename T>
SparseMatrixCSR<const T> SparseMatrixCSR<T>::to_const() const
{
    return const_type(*this);
}

template <typename T>
KOKKOS_FUNCTION int SparseMatrixCSR<T>::rows() const
{
    assert(this->rows_ >= 0);
    return this->rows_;
}
template <typename T>
KOKKOS_FUNCTION int SparseMatrixCSR<T>::columns() const
{
    assert(this->columns_ >= 0);
    return this->columns_;
}
template <typename T>
KOKKOS_FUNCTION int SparseMatrixCSR<T>::non_zero_size() const
{
    assert(this->nnz_ >= 0);
    assert(static_cast<size_t>(this->nnz_) <= static_cast<size_t>(this->rows_) * static_cast<size_t>(this->columns_));
    return this->nnz_;
}

template <typename T>
KOKKOS_FUNCTION int SparseMatrixCSR<T>::row_nz_size(int row) const
{
    assert(row >= 0 && row < rows_);
    return row_start_indices_(row + 1) - row_start_indices_(row);
}

template <typename T>
KOKKOS_FUNCTION SparseMatrixCSR<T>::int_element_type& SparseMatrixCSR<T>::row_nz_index(int row, int nz_index) const
{
    assert(row >= 0 && row < rows_);
    assert(nz_index >= 0 && nz_index < row_nz_size(row));
    return column_indices_(row_start_indices_(row) + nz_index);
}

template <typename T>
KOKKOS_FUNCTION T& SparseMatrixCSR<T>::row_nz_entry(int row, int nz_index) const
{
    assert(row >= 0 && row < rows_);
    assert(nz_index >= 0 && nz_index < row_nz_size(row));
    return values_(row_start_indices_(row) + nz_index);
}

template <typename T>
KOKKOS_FUNCTION T* SparseMatrixCSR<T>::values_data() const
{
    return values_.data();
}

template <typename T>
KOKKOS_FUNCTION SparseMatrixCSR<T>::int_element_type* SparseMatrixCSR<T>::column_indices_data() const
{
    return column_indices_.data();
}

template <typename T>
KOKKOS_FUNCTION SparseMatrixCSR<T>::int_element_type* SparseMatrixCSR<T>::row_start_indices_data() const
{
    return row_start_indices_.data();
}
} // namespace gmgpolar
