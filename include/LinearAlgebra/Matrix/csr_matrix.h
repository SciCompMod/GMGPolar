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

namespace sparse_csr_helpers {
		template<class MemorySpace>
void row_start_indices_from_nz_per_row(int& nnz, const Vector<int, MemorySpace>& nz_per_row, int rows, const Vector<int, MemorySpace>& row_start_indices)
{
	using ExecSpace = std::conditional_t<std::is_same_v<MemorySpace, Kokkos::HostSpace>, Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultExecutionSpace>;
    nnz = 0;
	    Kokkos::parallel_reduce(
        "compute nz_per_row", Kokkos::RangePolicy<ExecSpace>(0, 1),
        KOKKOS_LAMBDA(const std::size_t i, int& local_nnz) {
    for (int i = 0; i < rows; i++) {
        row_start_indices(i) = local_nnz;
        local_nnz += nz_per_row(i);
    }
    row_start_indices(rows) = local_nnz;
        },
        nnz);
	Kokkos::fence();

}
}

template <typename T, class MemorySpace = Kokkos::HostSpace>
class SparseMatrixCSR
{
public:
    using triplet_type = std::tuple<int, int, T>;

    // Default constructor.
    SparseMatrixCSR();

    // Copy constructor — A shallow copy is intentional here
    KOKKOS_DEFAULTED_FUNCTION SparseMatrixCSR(const SparseMatrixCSR& other) = default;

    // Move constructor — takes ownership instead of sharing.
    KOKKOS_FUNCTION SparseMatrixCSR(SparseMatrixCSR&& other) noexcept;

	SparseMatrixCSR(int rows, int columns, Vector<int, MemorySpace> nz_per_row);
    explicit SparseMatrixCSR(int rows, int columns, std::function<int(int)> nz_per_row);
    explicit SparseMatrixCSR(int rows, int columns, const std::vector<triplet_type>& entries);
    explicit SparseMatrixCSR(int rows, int columns, const std::vector<T>& values,
                             const std::vector<int>& column_indices, const std::vector<int>& row_start_indices);

    // Copy assignment is deleted.
    // Copying is only allowed through the copy constructor or copy() below,
    // making ownership transfer explicit and preventing accidental shallow
    // copies that silently share mutable storage.
    SparseMatrixCSR& operator=(const SparseMatrixCSR& other) = delete;

    // Move assignment — transfers ownership.
    SparseMatrixCSR& operator=(SparseMatrixCSR&& other) noexcept = default;

    SparseMatrixCSR copy() const;
    template<class TargetMemorySpace>
    std::conditional_t<std::is_same_v<MemorySpace, TargetMemorySpace>, const SparseMatrixCSR&, SparseMatrixCSR<T, TargetMemorySpace>> mirror_view_and_copy() const;

    KOKKOS_FUNCTION int rows() const;
    KOKKOS_FUNCTION int columns() const;
    KOKKOS_FUNCTION int non_zero_size() const;

    KOKKOS_FUNCTION int row_nz_size(int row) const;

    KOKKOS_FUNCTION const int& row_nz_index(int row, int nz_index) const;
    KOKKOS_FUNCTION void set_row_nz_index(int row, int nz_index, int row_nz_index) const;
    KOKKOS_FUNCTION void increase_row_nz_index(int row, int nz_index, int row_nz_index) const;

    KOKKOS_FUNCTION const T& row_nz_entry(int row, int nz_index) const;
    KOKKOS_FUNCTION void set_row_nz_entry(int row, int nz_index, T row_nz_entry) const;
    KOKKOS_FUNCTION void increase_row_nz_entry(int row, int nz_index, T row_nz_entry) const;

    KOKKOS_FUNCTION T* values_data() const;
    KOKKOS_FUNCTION int* column_indices_data() const;
    KOKKOS_FUNCTION int* row_start_indices_data() const;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const SparseMatrixCSR<U, MemorySpace>& matrix);

    template <class, class>
    friend class SparseMatrixCSR;

private:
    int rows_;
    int columns_;
    int nnz_;
    AllocatableVector<T, MemorySpace> values_;
    AllocatableVector<int, MemorySpace> column_indices_;
    AllocatableVector<int, MemorySpace> row_start_indices_;

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

template <typename U, class MemorySpace>
std::ostream& operator<<(std::ostream& stream, const SparseMatrixCSR<U, MemorySpace>& matrix)
{
    // Create host mirrors and deep_copy from device if necessary.
    // If the View is already on host, deep_copy is a no-op.
    auto h_values            = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, matrix.values_);
    auto h_column_indices    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, matrix.column_indices_);
    auto h_row_start_indices = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, matrix.row_start_indices_);

    stream << "SparseMatrixCSR: " << matrix.rows_ << " x " << matrix.columns_ << "\n";
    stream << "Number of non-zeros (nnz): " << matrix.nnz_ << "\n";
    stream << "Non-zero elements (row, column, value):\n";
    for (int row = 0; row < matrix.rows_; ++row) {
        for (int nnz = h_row_start_indices(row); nnz < h_row_start_indices(row + 1); ++nnz) {
            stream << "(" << row << ", " << h_column_indices(nnz) << ", " << h_values(nnz) << ")\n";
        }
    }
    return stream;
}

// default construction
template <typename T, class MemorySpace>
SparseMatrixCSR<T, MemorySpace>::SparseMatrixCSR()
    : rows_(0)
    , columns_(0)
    , nnz_(0)
{
}

// copy construction
template <typename T, class MemorySpace>
SparseMatrixCSR<T, MemorySpace> SparseMatrixCSR<T, MemorySpace>::copy() const
{
    SparseMatrixCSR<T, MemorySpace> new_copy;
    new_copy.rows_              = rows_;
    new_copy.columns_           = columns_;
    new_copy.nnz_               = nnz_;
    new_copy.values_            = AllocatableVector<T, MemorySpace>("CSR values", nnz_);
    new_copy.column_indices_    = AllocatableVector<int, MemorySpace>("CSR column indices", nnz_);
    new_copy.row_start_indices_ = AllocatableVector<int, MemorySpace>("CSR row start indices", rows_ + 1);
    Kokkos::deep_copy(new_copy.values_, values_);
    Kokkos::deep_copy(new_copy.column_indices_, column_indices_);
    Kokkos::deep_copy(new_copy.row_start_indices_, row_start_indices_);

    return new_copy;
}

template <typename T, class MemorySpace>
template <class TargetMemorySpace>
std::conditional_t<std::is_same_v<MemorySpace, TargetMemorySpace>, const SparseMatrixCSR<T, MemorySpace>&, SparseMatrixCSR<T, TargetMemorySpace>> SparseMatrixCSR<T, MemorySpace>::mirror_view_and_copy() const
{
if constexpr(std::is_same_v<MemorySpace, TargetMemorySpace>) {
		return *this;
} else {
    SparseMatrixCSR<T, TargetMemorySpace> new_copy;
    new_copy.rows_              = rows_;
    new_copy.columns_           = columns_;
    new_copy.nnz_               = nnz_;
    new_copy.values_            = Kokkos::create_mirror_view_and_copy(TargetMemorySpace(), values_);
    new_copy.column_indices_    = Kokkos::create_mirror_view_and_copy(TargetMemorySpace(), column_indices_);
    new_copy.row_start_indices_ = Kokkos::create_mirror_view_and_copy(TargetMemorySpace(), row_start_indices_);
    return new_copy;
}
}

// move construction
template <typename T, class MemorySpace>
KOKKOS_FUNCTION SparseMatrixCSR<T, MemorySpace>::SparseMatrixCSR(SparseMatrixCSR&& other) noexcept
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

template <typename T, class MemorySpace>
SparseMatrixCSR<T, MemorySpace>::SparseMatrixCSR(int rows, int columns, std::function<int(int)> nz_per_row)
    : rows_(rows)
    , columns_(columns)
    , row_start_indices_("CSR row start indices", rows_ + 1)
{
    assert(rows >= 0);
    assert(columns >= 0);

    auto h_row_start_indices = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, row_start_indices_);

    nnz_ = 0;
    for (int i = 0; i < rows; i++) {
        h_row_start_indices(i) = nnz_;
        nnz_ += nz_per_row(i);
    }
    h_row_start_indices(rows) = nnz_;
    values_                  = Vector<T, MemorySpace>("CSR values", nnz_);
    column_indices_          = Vector<int, MemorySpace>("CSR column indices", nnz_);

    assign(values_, T(0));
    assign(column_indices_, 0);
    Kokkos::deep_copy(row_start_indices_, h_row_start_indices);
}

template <typename T, class MemorySpace>
SparseMatrixCSR<T, MemorySpace>::SparseMatrixCSR(int rows, int columns, Vector<int, MemorySpace> nz_per_row)
    : rows_(rows)
    , columns_(columns)
    , row_start_indices_("CSR row start indices", rows + 1)
{
    assert(rows >= 0);
    assert(columns >= 0);

	sparse_csr_helpers::row_start_indices_from_nz_per_row(nnz_, nz_per_row, rows, row_start_indices_);
    values_                  = Vector<T, MemorySpace>("CSR values", nnz_);
    column_indices_          = Vector<int, MemorySpace>("CSR column indices", nnz_);

    assign(values_, T(0));
    assign(column_indices_, 0);
}

template <typename T, class MemorySpace>
SparseMatrixCSR<T, MemorySpace>::SparseMatrixCSR(int rows, int columns, const std::vector<triplet_type>& entries)
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

    auto h_values            = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, values_);
    auto h_column_indices    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, column_indices_);
    auto h_row_start_indices = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, row_start_indices_);

    // fill values and column indexes
    for (int i = 0; i < nnz_; i++) {
        assert(0 <= std::get<0>(entries[i]) && std::get<0>(entries[i]) < rows);
        h_values(i)         = std::get<2>(entries[i]);
        h_column_indices(i) = std::get<1>(entries[i]);
        assert(0 <= h_column_indices(i) && h_column_indices(i) < columns);
    }
    //fill row indexes
    int count             = 0;
    h_row_start_indices(0) = 0;
    for (int r = 0; r < rows; r++) {
        while (count < nnz_ && std::get<0>(entries[count]) == r)
            count++;
        h_row_start_indices(r + 1) = count;
    }
    assert(h_row_start_indices(rows) == nnz_);

    Kokkos::deep_copy(values_, h_values);
    Kokkos::deep_copy(column_indices_, h_column_indices);
    Kokkos::deep_copy(row_start_indices_, h_row_start_indices);
}

template <typename T, class MemorySpace>
SparseMatrixCSR<T, MemorySpace>::SparseMatrixCSR(int rows, int columns, const std::vector<T>& values,
                                                 const std::vector<int>& column_indices,
                                                 const std::vector<int>& row_start_indices)
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

    auto h_values            = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, values_);
    auto h_column_indices    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, column_indices_);
    auto h_row_start_indices = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, row_start_indices_);

    // Copy data to internal storage
    std::copy(values.begin(), values.end(), h_values.data());
    std::copy(column_indices.begin(), column_indices.end(), h_column_indices.data());
    std::copy(row_start_indices.begin(), row_start_indices.end(), h_row_start_indices.data());

    Kokkos::deep_copy(values_, h_values);
    Kokkos::deep_copy(column_indices_, h_column_indices);
    Kokkos::deep_copy(row_start_indices_, h_row_start_indices);
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION int SparseMatrixCSR<T, MemorySpace>::rows() const
{
    KOKKOS_ASSERT(this->rows_ >= 0);
    return this->rows_;
}
template <typename T, class MemorySpace>
KOKKOS_FUNCTION int SparseMatrixCSR<T, MemorySpace>::columns() const
{
    KOKKOS_ASSERT(this->columns_ >= 0);
    return this->columns_;
}
template <typename T, class MemorySpace>
KOKKOS_FUNCTION int SparseMatrixCSR<T, MemorySpace>::non_zero_size() const
{
    KOKKOS_ASSERT(this->nnz_ >= 0);
    KOKKOS_ASSERT(static_cast<size_t>(this->nnz_) <=
                  static_cast<size_t>(this->rows_) * static_cast<size_t>(this->columns_));
    return this->nnz_;
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION int SparseMatrixCSR<T, MemorySpace>::row_nz_size(int row) const
{
    KOKKOS_ASSERT(row >= 0 && row < rows_);
    return row_start_indices_(row + 1) - row_start_indices_(row);
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION const int& SparseMatrixCSR<T, MemorySpace>::row_nz_index(int row, int nz_index) const
{
    KOKKOS_ASSERT(row >= 0 && row < rows_);
    KOKKOS_ASSERT(nz_index >= 0 && nz_index < row_nz_size(row));
    return column_indices_(row_start_indices_(row) + nz_index);
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION void SparseMatrixCSR<T, MemorySpace>::set_row_nz_index(int row, int nz_index, int row_nz_index) const
{
    KOKKOS_ASSERT(row >= 0 && row < rows_);
    KOKKOS_ASSERT(nz_index >= 0 && nz_index < row_nz_size(row));
    column_indices_(row_start_indices_(row) + nz_index) = row_nz_index;
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION void SparseMatrixCSR<T, MemorySpace>::increase_row_nz_index(int row, int nz_index,
                                                                            int row_nz_index) const
{
    KOKKOS_ASSERT(row >= 0 && row < rows_);
    KOKKOS_ASSERT(nz_index >= 0 && nz_index < row_nz_size(row));
    column_indices_(row_start_indices_(row) + nz_index) += row_nz_index;
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION const T& SparseMatrixCSR<T, MemorySpace>::row_nz_entry(int row, int nz_index) const
{
    KOKKOS_ASSERT(row >= 0 && row < rows_);
    KOKKOS_ASSERT(nz_index >= 0 && nz_index < row_nz_size(row));
    return values_(row_start_indices_(row) + nz_index);
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION void SparseMatrixCSR<T, MemorySpace>::set_row_nz_entry(int row, int nz_index, T row_nz_entry) const
{
    KOKKOS_ASSERT(row >= 0 && row < rows_);
    KOKKOS_ASSERT(nz_index >= 0 && nz_index < row_nz_size(row));
    values_(row_start_indices_(row) + nz_index) = row_nz_entry;
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION void SparseMatrixCSR<T, MemorySpace>::increase_row_nz_entry(int row, int nz_index, T row_nz_entry) const
{
    KOKKOS_ASSERT(row >= 0 && row < rows_);
    KOKKOS_ASSERT(nz_index >= 0 && nz_index < row_nz_size(row));
    values_(row_start_indices_(row) + nz_index) += row_nz_entry;
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION T* SparseMatrixCSR<T, MemorySpace>::values_data() const
{
    return values_.data();
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION int* SparseMatrixCSR<T, MemorySpace>::column_indices_data() const
{
    return column_indices_.data();
}

template <typename T, class MemorySpace>
KOKKOS_FUNCTION int* SparseMatrixCSR<T, MemorySpace>::row_start_indices_data() const
{
    return row_start_indices_.data();
}
} // namespace gmgpolar
