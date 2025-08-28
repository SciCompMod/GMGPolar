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
#include <iomanip>

#include "csr_matrix.h"
#include "vector.h"

/**
 * @brief Sparse LU decomposition solver for symmetric positive definite matrices.
 *
 * This solver performs a sparse LU factorization with static pivoting and
 * diagonal perturbation to ensure numerical stability. It is slower than MUMPS,
 * but useful when an in-house LU implementation is desired.
 *
 * @tparam T Numeric type (e.g. double, float).
 */
template <typename T>
class SparseLUSolver
{
public:
    /**
     * @brief Construct an empty solver with given tolerances.
     *
     * @param tolerance_abs Minimum allowed diagonal magnitude.
     *        Any diagonal entry smaller than this will be perturbed.
     *        Default: 1e-12.
     * @param tolerance_rel Relative tolerance with respect to the largest
     *        entry in a row. Ensures diagonal is not too small compared
     *        to other row entries. Default: 1e-8.
     */
    explicit SparseLUSolver(T tolerance_abs = static_cast<T>(1e-12), T tolerance_rel = static_cast<T>(1e-8));

    /**
     * @brief Construct a solver from a sparse matrix and factorize it.
     *
     * @param A Square sparse matrix in CSR format to factorize.
     * @param tolerance_abs Minimum allowed diagonal magnitude.
     *        Any diagonal entry smaller than this will be perturbed.
     *        Default: 1e-12.
     * @param tolerance_rel Relative tolerance with respect to the largest
     *        entry in a row. Ensures diagonal is not too small compared
     *        to other row entries. Default: 1e-8.
     */
    explicit SparseLUSolver(const SparseMatrixCSR<T>& A, T tolerance_abs = static_cast<T>(1e-12),
                            T tolerance_rel = static_cast<T>(1e-8));

    SparseLUSolver(const SparseLUSolver& other);
    SparseLUSolver(SparseLUSolver&& other) noexcept;

    SparseLUSolver& operator=(const SparseLUSolver& other);
    SparseLUSolver& operator=(SparseLUSolver&& other) noexcept;

    /**
     * @brief Solve a linear system in place.
     *
     * This method overwrites the input vector with the solution `x`
     * to the system `Ax = b`, where `A` was the matrix provided at factorization.
     *
     * @param b Right-hand side vector (modified in place to contain the solution).
     */
    void solveInPlace(Vector<T>& b) const;
    void solveInPlace(T* b) const;

private:
    std::vector<T> L_values, U_values;
    std::vector<int> L_col_idx, U_col_idx;
    std::vector<int> L_row_ptr, U_row_ptr;
    bool factorized_ = false;

    T tolerance_abs_; // minimum allowed diagonal
    T tolerance_rel_; // relative to the max in the row

    void factorize(const SparseMatrixCSR<T>& A);

    void factorizeAccumulateSorted(const SparseMatrixCSR<T>& A);
    void factorizeWithHashing(const SparseMatrixCSR<T>& A);
};

// default construction
template <typename T>
SparseLUSolver<T>::SparseLUSolver(T tolerance_abs, T tolerance_rel)
    : factorized_(false)
    , tolerance_abs_(tolerance_abs)
    , tolerance_rel_(tolerance_rel)
{
}

// constructor
template <typename T>
SparseLUSolver<T>::SparseLUSolver(const SparseMatrixCSR<T>& A, T tolerance_abs, T tolerance_rel)
    : tolerance_abs_(tolerance_abs)
    , tolerance_rel_(tolerance_rel)
{
    assert(A.rows() == A.columns());
    if (!factorized_) {
        factorize(A);
    }
}

// copy construction
template <typename T>
SparseLUSolver<T>::SparseLUSolver(const SparseLUSolver& other)
    : L_values(other.L_values)
    , U_values(other.U_values)
    , L_col_idx(other.L_col_idx)
    , U_col_idx(other.U_col_idx)
    , L_row_ptr(other.L_row_ptr)
    , U_row_ptr(other.U_row_ptr)
    , factorized_(other.factorized_)
{
}

// copy assignment
template <typename T>
SparseLUSolver<T>& SparseLUSolver<T>::operator=(const SparseLUSolver& other)
{
    if (this == &other) {
        return *this; // Handle self-assignment
    }

    L_values    = other.L_values;
    U_values    = other.U_values;
    L_col_idx   = other.L_col_idx;
    U_col_idx   = other.U_col_idx;
    L_row_ptr   = other.L_row_ptr;
    U_row_ptr   = other.U_row_ptr;
    factorized_ = other.factorized_;

    return *this;
}

// move construction
template <typename T>
SparseLUSolver<T>::SparseLUSolver(SparseLUSolver&& other) noexcept
    : L_values(std::move(other.L_values))
    , U_values(std::move(other.U_values))
    , L_col_idx(std::move(other.L_col_idx))
    , U_col_idx(std::move(other.U_col_idx))
    , L_row_ptr(std::move(other.L_row_ptr))
    , U_row_ptr(std::move(other.U_row_ptr))
    , factorized_(other.factorized_)
{
    other.factorized_ = false;
}

// move assignment
template <typename T>
SparseLUSolver<T>& SparseLUSolver<T>::operator=(SparseLUSolver&& other) noexcept
{
    if (this == &other) {
        return *this; // Handle self-assignment
    }

    L_values    = std::move(other.L_values);
    U_values    = std::move(other.U_values);
    L_col_idx   = std::move(other.L_col_idx);
    U_col_idx   = std::move(other.U_col_idx);
    L_row_ptr   = std::move(other.L_row_ptr);
    U_row_ptr   = std::move(other.U_row_ptr);
    factorized_ = other.factorized_;

    other.factorized_ = false;

    return *this;
}

template <typename T>
void SparseLUSolver<T>::factorize(const SparseMatrixCSR<T>& A)
{
    factorizeWithHashing(A);
}

template <typename T>
void SparseLUSolver<T>::factorizeWithHashing(const SparseMatrixCSR<T>& A)
{
    const int n = A.rows();
    L_values.clear();
    U_values.clear();
    L_col_idx.clear();
    U_col_idx.clear();
    L_row_ptr.clear();
    U_row_ptr.clear();

    L_row_ptr.resize(n + 1, 0);
    U_row_ptr.resize(n + 1, 0);

    // Temporary structures to store computed rows
    std::vector<std::unordered_map<int, T>> L_map(n);
    std::vector<std::unordered_map<int, T>> U_map(n);

    for (int i = 0; i < n; i++) {
        std::unordered_map<int, T> row_values;

        // Load nonzero elements of row i
        for (int idx = 0; idx < A.row_nz_size(i); idx++) {
            int j         = A.row_nz_index(i, idx);
            row_values[j] = A.row_nz_entry(i, idx);
        }

        for (int j = 0; j < i; j++) {
            auto it = row_values.find(j);
            if (it == row_values.end())
                continue;

            // Compute L(i, j) using the diagonal U(j, j)
            it->second /= U_map[j][j];
            L_map[i][j] = it->second;

            // Update the remaining values in row i using the jth row of U.
            for (const auto& [k, U_val] : U_map[j]) {
                if (k > j) {
                    row_values[k] -= it->second * U_val;
                }
            }
        }

        for (const auto& [j, val] : row_values) {
            if (j < i)
                L_map[i][j] = val;
            else
                U_map[i][j] = val;
        }

        // Static pivoting / diagonal perturbation
        T& diag   = U_map[i][i]; // reference to the diagonal
        T max_val = 0;
        for (const auto& [col, val] : U_map[i]) {
            max_val = std::max(max_val, std::abs(val));
        }

        T threshold = std::max(tolerance_abs_, tolerance_rel_ * max_val);
        if (std::abs(diag) < threshold) {
            diag = std::copysign(threshold, diag);
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
void SparseLUSolver<T>::solveInPlace(T* b) const
{
    assert(factorized_);
    const int n = L_row_ptr.size() - 1; // n is the number of rows in the matrix

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
                diag = U_values[idx]; // Store diagonal value
            }
            else {
                b[i] -= U_values[idx] * b[col];
            }
        }
        assert(std::abs(diag) > tolerance_abs_);
        b[i] /= diag;
    }
}

template <typename T>
void SparseLUSolver<T>::solveInPlace(Vector<T>& b) const
{
    assert(b.size() == static_cast<int>(L_row_ptr.size()) - 1);
    solveInPlace(b.begin());
}
