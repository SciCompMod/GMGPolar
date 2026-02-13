#pragma once

#include <algorithm>
#include <queue>
#include <vector>
#include <cassert>
#include <cmath>
#include <stack>

#include "csr_matrix.h"
#include "vector.h"

/**
 * @brief Sparse LU decomposition solver for symmetric positive definite matrices.
 *
 * This solver performs a sparse LU factorization with static pivoting and
 * diagonal perturbation to ensure numerical stability. 
 * It is slightly slower than MUMPS, but useful when an in-house LU implementation is desired.
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

    /**
     * @brief Solve a linear system in place.
     *
     * This method overwrites the input vector with the solution `x`
     * to the system `Ax = b`, where `A` was the matrix provided at factorization.
     *
     * @param b Right-hand side vector (modified in place to contain the solution).
     */
    void solveInPlace(Vector<T> b) const;
    void solveInPlace(T* b) const;

private:
    // LU decomposition data structures
    std::vector<T> L_values, U_values; // Non-zero values for L and U
    std::vector<int> L_col_idx, U_col_idx; // Column indices for L and U
    std::vector<int> L_row_ptr, U_row_ptr; // Row pointers for L and U
    std::vector<T> U_diag; // Diagonal elements of U
    bool factorized_; // Factorization status flag
    T tolerance_abs_; // minimum allowed diagonal
    T tolerance_rel_; // relative to the max in the row

    // Core methods
    void factorize(const SparseMatrixCSR<T>& A);

    // Factorization components
    void symbolicFactorization(const SparseMatrixCSR<T>& A, std::vector<std::vector<int>>& L_pattern,
                               std::vector<std::vector<int>>& U_pattern) const;
    void numericFactorization(const SparseMatrixCSR<T>& A, const std::vector<std::vector<int>>& L_pattern,
                              const std::vector<std::vector<int>>& U_pattern);
};

// Default constructor
template <typename T>
SparseLUSolver<T>::SparseLUSolver(T tolerance_abs, T tolerance_rel)
    : factorized_(false)
    , tolerance_abs_(tolerance_abs)
    , tolerance_rel_(tolerance_rel)
{
}

/**
 * Constructs LU solver with RCM reordering and matrix factorization
 * @param A - Input matrix (must be square)
 */
template <typename T>
SparseLUSolver<T>::SparseLUSolver(const SparseMatrixCSR<T>& A, T tolerance_abs, T tolerance_rel)
    : factorized_(false)
    , tolerance_abs_(tolerance_abs)
    , tolerance_rel_(tolerance_rel)
{
    assert(A.rows() == A.columns());
    factorize(A);
    factorized_ = true;
}

/**
 * Solves Ax = b for Vector<T> type
 * @param b - Right-hand side vector (overwritten with solution)
 */
template <typename T>
void SparseLUSolver<T>::solveInPlace(Vector<T> b) const
{
    solveInPlace(b.data());
}

/**
 * Solves Ax = b for raw pointer
 * @param b - Right-hand side vector (overwritten with solution)
 */
template <typename T>
void SparseLUSolver<T>::solveInPlace(T* b) const
{
    const int n = L_row_ptr.size() - 1;

    // Forward substitution: L * y = b
    for (int i = 0; i < n; i++) {
        for (int idx = L_row_ptr[i]; idx < L_row_ptr[i + 1]; idx++) {
            b[i] -= L_values[idx] * b[L_col_idx[idx]];
        }
    }

    // Backward substitution: U * x = y
    for (int i = n - 1; i >= 0; i--) {
        for (int idx = U_row_ptr[i]; idx < U_row_ptr[i + 1]; idx++) {
            const int col = U_col_idx[idx];
            if (col != i) { // Skip diagonal (handled separately)
                b[i] -= U_values[idx] * b[col];
            }
        }
        b[i] /= U_diag[i]; // Divide by diagonal element
    }
}

/**
 * Main factorization driver
 * @param A - Permuted matrix to factorize
 */
template <typename T>
void SparseLUSolver<T>::factorize(const SparseMatrixCSR<T>& A)
{
    std::vector<std::vector<int>> L_pattern, U_pattern;
    symbolicFactorization(A, L_pattern, U_pattern);
    numericFactorization(A, L_pattern, U_pattern);
    factorized_ = true;
}

/**
 * Symbolic factorization to determine L/U sparsity patterns
 * @param A - Input matrix
 * @param L_pattern - Output pattern for L
 * @param U_pattern - Output pattern for U
 */
template <typename T>
void SparseLUSolver<T>::symbolicFactorization(const SparseMatrixCSR<T>& A, std::vector<std::vector<int>>& L_pattern,
                                              std::vector<std::vector<int>>& U_pattern) const
{
    const int n = A.rows();
    L_pattern.resize(n);
    U_pattern.resize(n);

    std::vector<int> marker(n, -1);
    std::vector<int> stk;

    for (int i = 0; i < n; ++i) {
        L_pattern[i].clear();
        U_pattern[i].clear();
        stk.clear();

        // Process original non-zeros in row i
        const int nnz = A.row_nz_size(i);
        for (int idx = 0; idx < nnz; ++idx) {
            int col = A.row_nz_index(i, idx);
            if (marker[col] != i) {
                marker[col] = i;
                if (col < i) {
                    L_pattern[i].push_back(col);
                    stk.push_back(col);
                }
                else {
                    U_pattern[i].push_back(col);
                }
            }
        }

        // Ensure diagonal is included in U
        if (marker[i] != i) {
            marker[i] = i;
            U_pattern[i].push_back(i);
        }

        // Process fill-ins
        while (!stk.empty()) {
            int j = stk.back();
            stk.pop_back();

            for (int u_col : U_pattern[j]) {
                if (marker[u_col] != i) {
                    marker[u_col] = i;
                    if (u_col < i) {
                        L_pattern[i].push_back(u_col);
                        stk.push_back(u_col);
                    }
                    else {
                        U_pattern[i].push_back(u_col);
                    }
                }
            }
        }

        // Sort patterns for efficient access
        std::sort(L_pattern[i].begin(), L_pattern[i].end());
        std::sort(U_pattern[i].begin(), U_pattern[i].end());
    }
}

/**
 * Numeric factorization using symbolic patterns
 * @param A - Input matrix
 * @param L_pattern - Symbolic pattern for L
 * @param U_pattern - Symbolic pattern for U
 */
template <typename T>
void SparseLUSolver<T>::numericFactorization(const SparseMatrixCSR<T>& A,
                                             const std::vector<std::vector<int>>& L_pattern,
                                             const std::vector<std::vector<int>>& U_pattern)
{
    const int n = A.rows();

    // Initialize storage structures
    L_row_ptr.resize(n + 1, 0);
    U_row_ptr.resize(n + 1, 0);
    U_diag.resize(n, 0);

    // Compute row pointers
    for (int i = 0; i < n; i++) {
        L_row_ptr[i + 1] = L_row_ptr[i] + L_pattern[i].size();
        U_row_ptr[i + 1] = U_row_ptr[i] + U_pattern[i].size();
    }

    // Allocate memory for values and indices
    L_values.resize(L_row_ptr[n]);
    L_col_idx.resize(L_row_ptr[n]);
    U_values.resize(U_row_ptr[n]);
    U_col_idx.resize(U_row_ptr[n]);

    // Find start of upper triangular part in U patterns
    std::vector<int> U_pattern_start_upper(n, 0);
    for (int j = 0; j < n; j++) {
        size_t pos = 0;
        while (pos < U_pattern[j].size() && U_pattern[j][pos] <= j)
            pos++;
        U_pattern_start_upper[j] = pos;
    }

    // Workspace for dense row computation
    std::vector<T> dense(n, 0);
    std::vector<int> marker(n, -1);
    std::stack<int> indices_used;

    for (int i = 0; i < n; i++) {
        // Initialize dense row with original matrix values
        for (int idx = 0; idx < A.row_nz_size(i); idx++) {
            int j    = A.row_nz_index(i, idx);
            T val    = A.row_nz_entry(i, idx);
            dense[j] = val;
            if (marker[j] != i) {
                marker[j] = i;
                indices_used.push(j);
            }
        }

        // Compute L elements
        int L_offset = L_row_ptr[i];
        for (const int j : L_pattern[i]) {
            const T Lij = dense[j] / U_diag[j];

            L_values[L_offset]  = Lij;
            L_col_idx[L_offset] = j;
            L_offset++;

            // Update dense row: dense -= Lij * U_row[j] (for columns k > j)
            const int U_update_start_offset = U_row_ptr[j] + U_pattern_start_upper[j];
            const int U_row_end_offset      = U_row_ptr[j + 1];
            const int update_len            = U_row_end_offset - U_update_start_offset;

            if (update_len <= 0) {
                continue;
            }
            const int* p_U_col = &U_col_idx[U_update_start_offset];
            const T* p_U_val   = &U_values[U_update_start_offset];

            // This inner loop is the most performance-critical part.
            for (int idx = 0; idx < update_len; ++idx) {
                const int k = p_U_col[idx];
                if (marker[k] != i) {
                    marker[k] = i;
                    dense[k]  = 0;
                    indices_used.push(k);
                }
                dense[k] -= Lij * p_U_val[idx];
            }
        }

        // Extract U elements and find diagonal
        T max_val           = 0;
        bool diagonal_found = false;
        int diag_offset     = -1;
        int U_offset        = U_row_ptr[i];
        for (int j : U_pattern[i]) {
            T val               = dense[j];
            U_values[U_offset]  = val;
            U_col_idx[U_offset] = j;
            if (j == i) {
                U_diag[i]      = val;
                diagonal_found = true;
                diag_offset    = U_offset;
            }
            // Track maximum value for pivoting
            T abs_val = std::abs(val);
            if (abs_val > max_val)
                max_val = abs_val;
            U_offset++;
        }

        // Static pivoting for weak diagonals
        T threshold_val = std::max(tolerance_abs_, tolerance_rel_ * max_val);
        if (std::abs(U_diag[i]) < threshold_val) {
            U_diag[i] = std::copysign(threshold_val, U_diag[i]);
            if (diag_offset != -1) {
                U_values[diag_offset] = U_diag[i];
            }
        }

        // Reset workspace for next row
        while (!indices_used.empty()) {
            int idx = indices_used.top();
            indices_used.pop();
            dense[idx]  = 0;
            marker[idx] = -1;
        }
    }
}
