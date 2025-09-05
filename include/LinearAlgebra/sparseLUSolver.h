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
 * We utilize Reverse Cuthill-McKee (RCM) reordering to minimize fill-in during the factorization process.
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
    void solveInPlace(Vector<T>& b) const;
    void solveInPlace(T* b) const;

private:
    // LU decomposition data structures
    std::vector<T> L_values, U_values; // Non-zero values for L and U
    std::vector<int> L_col_idx, U_col_idx; // Column indices for L and U
    std::vector<int> L_row_ptr, U_row_ptr; // Row pointers for L and U
    std::vector<int> perm; // Permutation vector (RCM ordering)
    std::vector<int> perm_inv; // Inverse permutation
    std::vector<T> U_diag; // Diagonal elements of U
    bool factorized_; // Factorization status flag
    T tolerance_abs_; // minimum allowed diagonal
    T tolerance_rel_; // relative to the max in the row

    // Core methods
    void factorize(const SparseMatrixCSR<T>& A);
    void solveInPlacePermuted(T* b) const;

    // Reordering and permutation utilities
    std::vector<int> computeRCM(const SparseMatrixCSR<T>& A) const;
    SparseMatrixCSR<T> permuteMatrix(const SparseMatrixCSR<T>& A, const std::vector<int>& perm,
                                     const std::vector<int>& perm_inv) const;

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

    // Compute RCM ordering
    perm = computeRCM(A);
    perm_inv.resize(perm.size());
    for (size_t i = 0; i < perm.size(); i++) {
        perm_inv[perm[i]] = i;
    }

    // Permute matrix according to RCM ordering
    SparseMatrixCSR<T> A_perm = permuteMatrix(A, perm, perm_inv);
    factorize(A_perm);
    factorized_ = true;
}

/**
 * Solves Ax = b for Vector<T> type
 * @param b - Right-hand side vector (overwritten with solution)
 */
template <typename T>
void SparseLUSolver<T>::solveInPlace(Vector<T>& b) const
{
    solveInPlace(b.begin());
}

/**
 * Solves Ax = b for raw pointer
 * @param b - Right-hand side vector (overwritten with solution)
 */
template <typename T>
void SparseLUSolver<T>::solveInPlace(T* b) const
{
    assert(factorized_);
    const int n = perm.size();
    if (n == 0)
        return;

    // Permute RHS: b_perm = P * b
    Vector<T> b_perm(n);
    for (int i = 0; i < n; i++) {
        b_perm[i] = b[perm[i]];
    }

    // Solve permuted system
    solveInPlacePermuted(b_perm.begin());

    // Unpermute solution: x = P^T * x_perm
    for (int i = 0; i < n; i++) {
        b[i] = b_perm[perm_inv[i]];
    }
}

/**
 * Performs forward/backward substitution on permuted system
 * @param b - Permuted right-hand side vector (overwritten with solution)
 */
template <typename T>
void SparseLUSolver<T>::solveInPlacePermuted(T* b) const
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
 * Computes Reverse Cuthill-McKee (RCM) ordering for bandwidth reduction
 * @param A - Input sparse matrix
 * @return Permutation vector
 */
template <typename T>
std::vector<int> SparseLUSolver<T>::computeRCM(const SparseMatrixCSR<T>& A) const
{
    const int n = A.rows();
    if (n == 0)
        return {};

    // Build symmetric adjacency list
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; i++) {
        for (int idx = 0; idx < A.row_nz_size(i); idx++) {
            int j = A.row_nz_index(i, idx);
            if (j == i)
                continue; // Skip diagonal
            adj[i].push_back(j);
            adj[j].push_back(i);
        }
    }

    // Remove duplicates and sort
    for (int i = 0; i < n; i++) {
        std::sort(adj[i].begin(), adj[i].end());
        auto last = std::unique(adj[i].begin(), adj[i].end());
        adj[i].resize(last - adj[i].begin());
    }

    // Compute degrees
    std::vector<int> deg(n);
    for (int i = 0; i < n; i++) {
        deg[i] = adj[i].size();
    }

    std::vector<bool> visited(n, false);
    std::vector<int> RCM_order;

    // Process disconnected components
    for (int start = 0; start < n; start++) {
        if (visited[start])
            continue;

        // Find connected component
        std::vector<int> comp_nodes;
        std::queue<int> q_comp;
        q_comp.push(start);
        visited[start] = true;
        while (!q_comp.empty()) {
            int u = q_comp.front();
            q_comp.pop();
            comp_nodes.push_back(u);
            for (int v : adj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    q_comp.push(v);
                }
            }
        }

        // Find min-degree node in component
        int comp_start = comp_nodes[0];
        int min_deg    = deg[comp_start];
        for (int node : comp_nodes) {
            visited[node] = false; // Unmark for BFS
            if (deg[node] < min_deg) {
                min_deg    = deg[node];
                comp_start = node;
            }
        }

        // BFS traversal
        std::queue<int> q;
        std::vector<int> comp_order;
        q.push(comp_start);
        visited[comp_start] = true;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            comp_order.push_back(u);

            // Collect unvisited neighbors
            std::vector<int> neighbors;
            for (int v : adj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    neighbors.push_back(v);
                }
            }

            // Sort neighbors by increasing degree
            std::sort(neighbors.begin(), neighbors.end(), [&](int a, int b) {
                return deg[a] < deg[b];
            });

            for (int v : neighbors) {
                q.push(v);
            }
        }

        // Reverse for RCM ordering and append
        std::reverse(comp_order.begin(), comp_order.end());
        RCM_order.insert(RCM_order.end(), comp_order.begin(), comp_order.end());
    }

    return RCM_order;
}

/**
 * Permutes matrix using RCM ordering (efficient, no sorting)
 * @param A - Original matrix
 * @param perm - Permutation vector
 * @param perm_inv - Inverse permutation
 * @return Permuted CSR matrix
 */
template <typename T>
SparseMatrixCSR<T> SparseLUSolver<T>::permuteMatrix(const SparseMatrixCSR<T>& A, const std::vector<int>& perm,
                                                    const std::vector<int>& perm_inv) const
{
    const int n = A.rows();

    // Compute number of nonzeros per permuted row
    std::vector<int> nz_per_row(n);
    for (int i_new = 0; i_new < n; ++i_new) {
        int i_old         = perm[i_new];
        nz_per_row[i_new] = A.row_nz_size(i_old);
    }

    // Construct permuted matrix with preallocated storage
    SparseMatrixCSR<T> A_perm(n, n, [&](int i) {
        return nz_per_row[i];
    });

    // Fill values and column indices
    for (int i_new = 0; i_new < n; ++i_new) {
        int i_old = perm[i_new];
        int nnz   = A.row_nz_size(i_old);
        for (int idx = 0; idx < nnz; ++idx) {
            int j_old = A.row_nz_index(i_old, idx);
            T val     = A.row_nz_entry(i_old, idx);
            int j_new = perm_inv[j_old];

            // Find the position in the underlying storage
            A_perm.row_nz_entry(i_new, idx) = val;
            A_perm.row_nz_index(i_new, idx) = j_new;
        }
    }

    return A_perm;
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
    L_pattern.clear();
    U_pattern.clear();
    L_pattern.resize(n);
    U_pattern.resize(n);

    // Marker array tracks visited columns per row
    std::vector<int> marker(n, -1);
    std::vector<int> row_marked_indices;
    std::vector<int> stk;

    for (int i = 0; i < n; ++i) {
        row_marked_indices.clear();
        stk.clear();

        // Process original non-zeros in row i
        const int nnz = A.row_nz_size(i);
        for (int idx = 0; idx < nnz; ++idx) {
            int col = A.row_nz_index(i, idx);
            if (marker[col] != i) {
                marker[col] = i;
                row_marked_indices.push_back(col);
                if (col < i) { // Lower triangular element
                    stk.push_back(col);
                }
            }
        }

        // Ensure diagonal is included
        if (marker[i] != i) {
            marker[i] = i;
            row_marked_indices.push_back(i);
        }

        // Process fill-in elements
        while (!stk.empty()) {
            int j = stk.back();
            stk.pop_back();

            // Process U[j] pattern
            for (int col2 : U_pattern[j]) {
                if (marker[col2] != i) {
                    marker[col2] = i;
                    row_marked_indices.push_back(col2);
                    if (col2 < i) { // Continue for lower indices
                        stk.push_back(col2);
                    }
                }
            }
        }

        // Split into L and U patterns
        auto& Li = L_pattern[i];
        auto& Ui = U_pattern[i];
        for (int col : row_marked_indices) {
            if (col < i) {
                Li.push_back(col);
            }
            else {
                Ui.push_back(col);
            }
        }

        // Sort patterns for efficient access
        std::sort(Li.begin(), Li.end());
        std::sort(Ui.begin(), Ui.end());
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
