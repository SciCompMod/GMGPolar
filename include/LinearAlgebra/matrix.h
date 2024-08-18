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
class SparseMatrix
{
public:
    using triplet_type = std::tuple<int, int, T>;

    SparseMatrix();
    SparseMatrix(const SparseMatrix& other);
    SparseMatrix(SparseMatrix&& other) noexcept;

    explicit SparseMatrix(int rows, int columns, int nnz);
    explicit SparseMatrix(int rows, int columns, const std::vector<triplet_type>& entries);

    SparseMatrix& operator=(const SparseMatrix& other);
    SparseMatrix& operator=(SparseMatrix&& other) noexcept;

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

    template<typename U>
    friend std::ostream& operator<<(std::ostream& stream, const SparseMatrix<U>& matrix){
        stream << "SparseMatrix: " << matrix.rows_ << " x " << matrix.columns_ << "\n";
        stream << "Number of non-zeros (nnz): " << matrix.nnz_ << "\n";
        if (matrix.is_symmetric_) {
            stream << "Matrix is symmetric.\n";
        }
        stream << "Non-zero elements (row, column, value):\n";
        for (int i = 0; i < matrix.nnz_; ++i) {
            stream << "(" << matrix.row_indices_[i] << ", "
                        << matrix.column_indices_[i] << ", "
                        << matrix.values_[i] << ")\n";
        }
        return stream;
    }
    
    void write_to_file(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open file");
        }
        file << "SparseMatrix: " << rows_ << " x " << columns_ << "\n";
        file << "Number of non-zeros (nnz): " << nnz_ << "\n";
        if (is_symmetric_) {
            file << "Matrix is symmetric.\n";
        }
        file << "Non-zero elements (row, column, value):\n";
        for (int i = 0; i < nnz_; ++i) {
            file << "(" << row_indices_[i] << ", " << column_indices_[i] << ", " << values_[i] << ")\n";
        }
        file.close();
    }

private:
    int rows_;
    int columns_;
    int nnz_;
    std::unique_ptr<int[]> row_indices_;
    std::unique_ptr<int[]> column_indices_;
    std::unique_ptr<T[]> values_;
    bool is_symmetric_ = false;
};

template<typename T>
void sort_entries(std::vector<std::tuple<int, int, T>>& entries)
{
    const auto compare = [](const auto entry1, const auto entry2)
    {
        const auto local_r1 = std::get<0>(entry1);
        const auto local_r2 = std::get<0>(entry2);
        if(local_r1 < local_r2)
        {
            return true;
        }
        else if (local_r1 == local_r2)
        {
            return std::get<1>(entry1) < std::get<1>(entry2);
        }
        return false;
    };
    std::sort(entries.begin(), entries.end(), compare);
}

// default construction
template<typename T>
SparseMatrix<T>::SparseMatrix() :
    rows_(0),
    columns_(0),
    nnz_(0),
    row_indices_(nullptr),
    column_indices_(nullptr),
    values_(nullptr),
    is_symmetric_(false)
{}

// copy construction
template<typename T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix& other) :
    rows_(other.rows_),
    columns_(other.columns_),
    nnz_(other.nnz_),
    row_indices_(std::make_unique<int[]>(nnz_)),
    column_indices_(std::make_unique<int[]>(nnz_)),
    values_(std::make_unique<T[]>(nnz_)),
    is_symmetric_(other.is_symmetric_)
{
    std::copy(other.row_indices_.get(), other.row_indices_.get() + nnz_, row_indices_.get());
    std::copy(other.column_indices_.get(), other.column_indices_.get() + nnz_, column_indices_.get());
    std::copy(other.values_.get(), other.values_.get() + nnz_, values_.get());
}

// copy assignment
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator=(const SparseMatrix& other){
    if (this == &other) {
        // Self-assignment, no work needed
        return *this;
    }
    // Only allocate new memory if the sizes are different
    if (nnz_ != other.nnz_) {
        row_indices_ = std::make_unique<int[]>(nnz_);
        column_indices_ = std::make_unique<int[]>(nnz_);
        values_ = std::make_unique<T[]>(nnz_);
    }
    // Copy the elements
    rows_ = other.rows_;
    columns_= other.columns_;
    nnz_ = other.nnz_;
    is_symmetric_ = other.is_symmetric_;
    std::copy(other.row_indices_.get(), other.row_indices_.get() + nnz_, row_indices_.get());
    std::copy(other.column_indices_.get(), other.column_indices_.get() + nnz_, column_indices_.get());
    std::copy(other.values_.get(), other.values_.get() + nnz_, values_.get());
    return *this;
}

// move construction
template<typename T>
SparseMatrix<T>::SparseMatrix(SparseMatrix&& other) noexcept :
    rows_(other.rows_),
    columns_(other.columns_),
    nnz_(other.nnz_),
    row_indices_(std::move(other.row_indices_)),
    column_indices_(std::move(other.column_indices_)),
    values_(std::move(other.values_)),
    is_symmetric_(other.is_symmetric_)
{
    other.nnz_ = 0;
    other.rows_ = 0;
    other.columns_ = 0;
    other.is_symmetric_ = false;
}

// move assignment
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator=(SparseMatrix&& other) noexcept{
    rows_ = other.rows_;
    columns_ = other.columns_;
    nnz_ = other.nnz_;
    row_indices_ = std::move(other.row_indices_);
    column_indices_ = std::move(other.column_indices_);
    values_ = std::move(other.values_);
    is_symmetric_ = other.is_symmetric_;
    other.nnz_ = 0;
    other.rows_ = 0;
    other.columns_ = 0;
    other.is_symmetric_ = false;
    return *this;
}

template<typename T>
SparseMatrix<T>::SparseMatrix(int rows, int columns, int nnz):
    rows_(rows),
    columns_(columns),
    nnz_(nnz),
    row_indices_(std::make_unique<int[]>(nnz)),
    column_indices_(std::make_unique<int[]>(nnz)),
    values_(std::make_unique<T[]>(nnz)),
    is_symmetric_(false)
{
    assert(rows >= 0);
    assert(columns >= 0);
    assert(nnz >= 0);
}

template<typename T>
SparseMatrix<T>::SparseMatrix(int rows, int columns, const std::vector<triplet_type>& entries):
    // entries: row_idx, col_idx, value
    rows_(rows),
    columns_(columns),
    nnz_(entries.size()),
    row_indices_(std::make_unique<int[]>(nnz_)),
    column_indices_(std::make_unique<int[]>(nnz_)),
    values_(std::make_unique<T[]>(nnz_)),
    is_symmetric_(false)
{
    assert(rows_ >= 0);
    assert(columns_ >= 0);
    assert(nnz_ >= 0);
    #pragma omp parallel for
    for (int i = 0; i < nnz_; i++){
        assert(0 <= std::get<0>(entries[i]) && std::get<0>(entries[i]) < rows_);
        assert(0 <= std::get<1>(entries[i]) && std::get<1>(entries[i]) < columns_);
        row_indices_[i] = std::get<0>(entries[i]);
        column_indices_[i] = std::get<1>(entries[i]);
        values_[i] = std::get<2>(entries[i]);
    }
}

template<typename T>
int SparseMatrix<T>::rows() const {
    assert(this->rows_ >= 0);
    return this->rows_;
}
template<typename T>
int SparseMatrix<T>::columns() const {
    assert(this->columns_ >= 0);
    return this->columns_;
}
template<typename T>
int SparseMatrix<T>::non_zero_size() const {
    assert(this->nnz_ >= 0);
    assert(static_cast<size_t>(this->nnz_) <= static_cast<size_t>(this->rows_) * static_cast<size_t>(this->columns_));
    return this->nnz_;
}

template<typename T>
int& SparseMatrix<T>::row_index(int nz_index) {
    assert(nz_index >= 0); assert(nz_index < this->nnz_);
    return this->row_indices_[nz_index];
}
template<typename T>
const int& SparseMatrix<T>::row_index(int nz_index) const {
    assert(nz_index >= 0); assert(nz_index < this->nnz_);
    return this->row_indices_[nz_index];
}

template<typename T>
int& SparseMatrix<T>::col_index(int nz_index) {
    assert(nz_index >= 0); assert(nz_index < this->nnz_);
    return this->column_indices_[nz_index];
}
template<typename T>
const int& SparseMatrix<T>::col_index(int nz_index) const {
    assert(nz_index >= 0); assert(nz_index < this->nnz_);
    return this->column_indices_[nz_index];
}

template<typename T>
T& SparseMatrix<T>::value(int nz_index) {
    assert(nz_index >= 0); assert(nz_index < this->nnz_);
    return this->values_[nz_index];
}
template<typename T>
const T& SparseMatrix<T>::value(int nz_index) const {
    assert(nz_index >= 0); assert(nz_index < this->nnz_);
    return this->values_[nz_index];
}

template<typename T>
bool SparseMatrix<T>::is_symmetric() const{
    return is_symmetric_;
}
template<typename T>
void SparseMatrix<T>::is_symmetric(bool value){
    is_symmetric_ = value;
}

template<typename T>
int* SparseMatrix<T>::row_indices_data() const {
    return row_indices_.get();
}

template<typename T>
int* SparseMatrix<T>::column_indices_data() const {
    return column_indices_.get();
}

template<typename T>
T* SparseMatrix<T>::values_data() const {
    return values_.get();
}