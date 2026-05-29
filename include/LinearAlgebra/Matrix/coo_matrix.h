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

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

namespace gmgpolar
{

template <typename T, class MemorySpace = DefaultMemorySpace>
class SparseMatrixCOO
{
public:
    using triplet_type = std::tuple<int, int, T>;

    SparseMatrixCOO();
    KOKKOS_DEFAULTED_FUNCTION SparseMatrixCOO(const SparseMatrixCOO& other) = default;
    SparseMatrixCOO(SparseMatrixCOO&& other) noexcept;

    explicit SparseMatrixCOO(int rows, int columns, int nnz);
    explicit SparseMatrixCOO(int rows, int columns, const std::vector<triplet_type>& entries);

    SparseMatrixCOO& operator=(const SparseMatrixCOO& other) = delete;
    SparseMatrixCOO& operator=(SparseMatrixCOO&& other) noexcept;

    SparseMatrixCOO copy() const;
    template <class TargetMemorySpace>
    std::conditional_t<std::is_same_v<MemorySpace, TargetMemorySpace>, const SparseMatrixCOO<T, TargetMemorySpace>&,
                       SparseMatrixCOO<T, TargetMemorySpace>>
    mirror_view_and_copy() const;

    int rows() const;
    int columns() const;
    int non_zero_size() const;

    KOKKOS_INLINE_FUNCTION const int& row_index(int nz_index) const;
    KOKKOS_INLINE_FUNCTION void set_row_index(int nz_index, int row_index) const;
    KOKKOS_INLINE_FUNCTION void increment_row_index(int nz_index) const;

    KOKKOS_INLINE_FUNCTION const int& col_index(int nz_index) const;
    KOKKOS_INLINE_FUNCTION void set_col_index(int nz_index, int col_index) const;
    KOKKOS_INLINE_FUNCTION void increment_col_index(int nz_index) const;

    KOKKOS_INLINE_FUNCTION const T& value(int nz_index) const;
    KOKKOS_INLINE_FUNCTION void set_value(int nz_index, T value) const;
    KOKKOS_INLINE_FUNCTION void increase_value(int nz_index, T value) const;

    bool is_symmetric() const;
    void is_symmetric(bool value);

    int* row_indices_data() const;
    int* column_indices_data() const;
    T* values_data() const;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const SparseMatrixCOO<U>& matrix);

    // SparseMatrixCOO is a friend to versions of itself on different memory spaces to simplify the implementation of the copy operator
    template <class, class>
    friend class SparseMatrixCOO;

    void write_to_file(const std::string& filename) const;

private:
    int rows_;
    int columns_;
    int nnz_;
    AllocatableVector<int, MemorySpace> row_indices_;
    AllocatableVector<int, MemorySpace> column_indices_;
    AllocatableVector<T, MemorySpace> values_;
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

template <typename T, class MemorySpace>
void SparseMatrixCOO<T, MemorySpace>::write_to_file(const std::string& filename) const
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

template <typename T, class MemorySpace>
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
template <typename T, class MemorySpace>
SparseMatrixCOO<T, MemorySpace>::SparseMatrixCOO()
    : rows_(0)
    , columns_(0)
    , nnz_(0)
    , is_symmetric_(false)
{
}

// copy construction
template <typename T, class MemorySpace>
SparseMatrixCOO<T, MemorySpace> SparseMatrixCOO<T, MemorySpace>::copy() const
{
    SparseMatrixCOO<T, MemorySpace> other;
    other.rows_           = rows_;
    other.columns_        = columns_;
    other.nnz_            = nnz_;
    other.row_indices_    = Vector<int, MemorySpace>("COO row indices", nnz_);
    other.column_indices_ = Vector<int, MemorySpace>("COO column indices", nnz_);
    other.values_         = Vector<T, MemorySpace>("COO values", nnz_);
    other.is_symmetric_   = is_symmetric_;

    Kokkos::deep_copy(other.row_indices_, row_indices_);
    Kokkos::deep_copy(other.column_indices_, column_indices_);
    Kokkos::deep_copy(other.values_, values_);

    return other;
}

template <typename T, class MemorySpace>
template <class TargetMemorySpace>
std::conditional_t<std::is_same_v<MemorySpace, TargetMemorySpace>, const SparseMatrixCOO<T, TargetMemorySpace>&,
                   SparseMatrixCOO<T, TargetMemorySpace>>
SparseMatrixCOO<T, MemorySpace>::mirror_view_and_copy() const
{
    if constexpr (std::is_same_v<MemorySpace, TargetMemorySpace>) {
        return *this;
    }
    else {
        SparseMatrixCOO<T, TargetMemorySpace> matrix(rows_, columns_, nnz_);
        matrix.row_indices_    = Kokkos::create_mirror_view_and_copy(TargetMemorySpace(), row_indices_);
        matrix.column_indices_ = Kokkos::create_mirror_view_and_copy(TargetMemorySpace(), column_indices_);
        matrix.values_         = Kokkos::create_mirror_view_and_copy(TargetMemorySpace(), values_);
        matrix.is_symmetric(is_symmetric_);
        return matrix;
    }
}

// move construction
template <typename T, class MemorySpace>
SparseMatrixCOO<T, MemorySpace>::SparseMatrixCOO(SparseMatrixCOO&& other) noexcept
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
template <typename T, class MemorySpace>
SparseMatrixCOO<T, MemorySpace>& SparseMatrixCOO<T, MemorySpace>::operator=(SparseMatrixCOO&& other) noexcept
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

template <typename T, class MemorySpace>
SparseMatrixCOO<T, MemorySpace>::SparseMatrixCOO(int rows, int columns, int nnz)
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

    assign(row_indices_, 0);
    assign(column_indices_, 0);
    assign(values_, T(0));
}

template <typename T, class MemorySpace>
SparseMatrixCOO<T, MemorySpace>::SparseMatrixCOO(int rows, int columns, const std::vector<triplet_type>& entries)
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
    auto h_row_indices    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, row_indices_);
    auto h_column_indices = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, column_indices_);
    auto h_values         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, values_);
    Kokkos::parallel_for("SparseMatrixCOO constructor", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, nnz_),
                         [&](const int i) {
                             assert(0 <= std::get<0>(entries[i]) && std::get<0>(entries[i]) < rows_);
                             assert(0 <= std::get<1>(entries[i]) && std::get<1>(entries[i]) < columns_);
                             h_row_indices(i)    = std::get<0>(entries[i]);
                             h_column_indices(i) = std::get<1>(entries[i]);
                             h_values(i)         = std::get<2>(entries[i]);
                         });
    Kokkos::deep_copy(row_indices_, h_row_indices);
    Kokkos::deep_copy(column_indices_, h_column_indices);
    Kokkos::deep_copy(values_, h_values);
}

template <typename T, class MemorySpace>
int SparseMatrixCOO<T, MemorySpace>::rows() const
{
    assert(this->rows_ >= 0);
    return this->rows_;
}
template <typename T, class MemorySpace>
int SparseMatrixCOO<T, MemorySpace>::columns() const
{
    assert(this->columns_ >= 0);
    return this->columns_;
}
template <typename T, class MemorySpace>
int SparseMatrixCOO<T, MemorySpace>::non_zero_size() const
{
    assert(this->nnz_ >= 0);
    assert(static_cast<size_t>(this->nnz_) <= static_cast<size_t>(this->rows_) * static_cast<size_t>(this->columns_));
    return this->nnz_;
}

template <typename T, class MemorySpace>
KOKKOS_INLINE_FUNCTION void SparseMatrixCOO<T, MemorySpace>::set_row_index(int nz_index, int row_index) const
{
    assert(nz_index >= 0);
    assert(nz_index < this->nnz_);
    row_indices_(nz_index) = row_index;
}
template <typename T, class MemorySpace>
KOKKOS_INLINE_FUNCTION const int& SparseMatrixCOO<T, MemorySpace>::row_index(int nz_index) const
{
    assert(nz_index >= 0);
    assert(nz_index < this->nnz_);
    return row_indices_(nz_index);
}
template <typename T, class MemorySpace>
KOKKOS_INLINE_FUNCTION void SparseMatrixCOO<T, MemorySpace>::increment_row_index(int nz_index) const
{
    assert(nz_index >= 0);
    assert(nz_index < this->nnz_);
    row_indices_(nz_index)++;
}

template <typename T, class MemorySpace>
KOKKOS_INLINE_FUNCTION void SparseMatrixCOO<T, MemorySpace>::set_col_index(int nz_index, int col_index) const
{
    assert(nz_index >= 0);
    assert(nz_index < nnz_);
    column_indices_(nz_index) = col_index;
}
template <typename T, class MemorySpace>
KOKKOS_INLINE_FUNCTION const int& SparseMatrixCOO<T, MemorySpace>::col_index(int nz_index) const
{
    assert(nz_index >= 0);
    assert(nz_index < nnz_);
    return column_indices_(nz_index);
}
template <typename T, class MemorySpace>
KOKKOS_INLINE_FUNCTION void SparseMatrixCOO<T, MemorySpace>::increment_col_index(int nz_index) const
{
    assert(nz_index >= 0);
    assert(nz_index < nnz_);
    column_indices_(nz_index)++;
}

template <typename T, class MemorySpace>
KOKKOS_INLINE_FUNCTION void SparseMatrixCOO<T, MemorySpace>::set_value(int nz_index, T value) const
{
    assert(nz_index >= 0);
    assert(nz_index < nnz_);
    values_(nz_index) = value;
}
template <typename T, class MemorySpace>
KOKKOS_INLINE_FUNCTION const T& SparseMatrixCOO<T, MemorySpace>::value(int nz_index) const
{
    assert(nz_index >= 0);
    assert(nz_index < nnz_);
    return values_(nz_index);
}
template <typename T, class MemorySpace>
KOKKOS_INLINE_FUNCTION void SparseMatrixCOO<T, MemorySpace>::increase_value(int nz_index, T value) const
{
    assert(nz_index >= 0);
    assert(nz_index < nnz_);
    values_(nz_index) += value;
}

template <typename T, class MemorySpace>
bool SparseMatrixCOO<T, MemorySpace>::is_symmetric() const
{
    return is_symmetric_;
}
template <typename T, class MemorySpace>
void SparseMatrixCOO<T, MemorySpace>::is_symmetric(bool value)
{
    is_symmetric_ = value;
}

template <typename T, class MemorySpace>
int* SparseMatrixCOO<T, MemorySpace>::row_indices_data() const
{
    return row_indices_.data();
}

template <typename T, class MemorySpace>
int* SparseMatrixCOO<T, MemorySpace>::column_indices_data() const
{
    return column_indices_.data();
}

template <typename T, class MemorySpace>
T* SparseMatrixCOO<T, MemorySpace>::values_data() const
{
    return values_.data();
}
} // namespace gmgpolar
