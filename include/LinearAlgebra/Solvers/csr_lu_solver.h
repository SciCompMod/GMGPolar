#pragma once

#include <algorithm>
#include <queue>
#include <vector>
#include <cassert>
#include <cmath>
#include <stack>

#include "../../LinearAlgebra/Matrix/csr_matrix.h"
#include "../../LinearAlgebra/Vector/vector.h"

namespace gmgpolar
{

namespace sparse_lu_helpers {

template <class MemorySpace>
static inline Vector<int, MemorySpace> build_perm_inv(Vector<int, MemorySpace> const& perm) {
	using ExecSpace = std::conditional_t<std::is_same_v<MemorySpace, Kokkos::HostSpace>, Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultExecutionSpace>;
    Vector<int, MemorySpace> perm_inv("perm_inv", perm.size());
    Kokkos::parallel_for(
        "Calculate perm_inv", Kokkos::RangePolicy<ExecSpace>(0, perm.size()),
        KOKKOS_LAMBDA(const int i) { perm_inv[perm[i]] = i; });
	Kokkos::fence();
    return perm_inv;
}

}

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

template <typename T, class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space>
class SparseLUSolver
{
	using ExecSpace = std::conditional_t<std::is_same_v<MemorySpace, Kokkos::HostSpace>, Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultExecutionSpace>;
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
    explicit SparseLUSolver(const SparseMatrixCSR<T, MemorySpace>& A, T tolerance_abs = static_cast<T>(1e-12),
                            T tolerance_rel = static_cast<T>(1e-8));

    /**
     * @brief Solve a linear system in place.
     *
     * This method overwrites the input vector with the solution `x`
     * to the system `Ax = b`, where `A` was the matrix provided at factorization.
     *
     * @param b Right-hand side vector (modified in place to contain the solution).
     */
    void solveInPlace(Vector<T, MemorySpace> b) const;

private:
    // LU decomposition data structures
    AllocatableVector<T, MemorySpace> L_values_, U_values_; // Non-zero values for L and U
    AllocatableVector<int, MemorySpace> L_col_idx_, U_col_idx_; // Column indices for L and U
    AllocatableVector<int, MemorySpace> L_row_ptr_, U_row_ptr_; // Row pointers for L and U
    AllocatableVector<int, MemorySpace> perm_; // Permutation vector (RCM ordering)
    AllocatableVector<int, MemorySpace> perm_inv_; // Inverse permutation
    AllocatableVector<T, MemorySpace> U_diag_; // Diagonal elements of U
    bool factorized_; // Factorization status flag
    T tolerance_abs_; // minimum allowed diagonal
    T tolerance_rel_; // relative to the max in the row

    // Core methods
    void factorize(const SparseMatrixCSR<T, MemorySpace>& A);
// Public due to cuda restrictions
public:
    void solveInPlacePermuted(const Vector<T, MemorySpace>& b) const;

private:
    // Reordering and permutation utilities
    Vector<int, MemorySpace> computeRCM(const SparseMatrixCSR<T, MemorySpace>& A) const;
// Public due to cuda restrictions
public:
    SparseMatrixCSR<T, MemorySpace> permuteMatrix(const SparseMatrixCSR<T, MemorySpace>& A,
                                                  const Vector<int, MemorySpace>& perm,
                                                  const Vector<int, MemorySpace>& perm_inv) const;

private:
    // Factorization components
    void symbolicFactorization(const SparseMatrixCSR<T, MemorySpace>& A, std::vector<std::vector<int>>& L_pattern,
                               std::vector<std::vector<int>>& U_pattern) const;
    void numericFactorization(const SparseMatrixCSR<T, MemorySpace>& A, const std::vector<std::vector<int>>& L_pattern,
                              const std::vector<std::vector<int>>& U_pattern);
};

// Default constructor
template <typename T, class MemorySpace>
SparseLUSolver<T, MemorySpace>::SparseLUSolver(T tolerance_abs, T tolerance_rel)
    : factorized_(false)
    , tolerance_abs_(tolerance_abs)
    , tolerance_rel_(tolerance_rel)
{
}

/**
 * Constructs LU solver with RCM reordering and matrix factorization
 * @param A - Input matrix (must be square)
 */
template <typename T, class MemorySpace>
SparseLUSolver<T, MemorySpace>::SparseLUSolver(const SparseMatrixCSR<T, MemorySpace>& A, T tolerance_abs,
                                               T tolerance_rel)
    : factorized_(false)
    , tolerance_abs_(tolerance_abs)
    , tolerance_rel_(tolerance_rel)
{
    assert(A.rows() == A.columns());

    // Compute RCM ordering
    perm_     = computeRCM(A);
    perm_inv_ = sparse_lu_helpers::build_perm_inv(perm_);

    // Permute matrix according to RCM ordering
    SparseMatrixCSR<T, MemorySpace> A_perm = permuteMatrix(A, perm_, perm_inv_);
    factorize(A_perm);
    factorized_ = true;
}

/**
 * Solves Ax = b for Vector<T> type
 * @param b - Right-hand side vector (overwritten with solution)
 */
template <typename T, class MemorySpace>
void SparseLUSolver<T, MemorySpace>::solveInPlace(Vector<T, MemorySpace> b) const
{
    assert(factorized_);
    const int n = perm_.size();
    if (n == 0)
        return;

	// Create local accessors to avoid copying the whole class to GPU
    Vector<int, MemorySpace> perm = perm_;
    Vector<int, MemorySpace> perm_inv = perm_inv_;

    // Permute RHS: b_perm = P * b
    Vector<T, MemorySpace> b_perm("b_perm", n);
    Kokkos::parallel_for(
        "b permute", Kokkos::RangePolicy<ExecSpace>(0, n), KOKKOS_LAMBDA(const int i) { b_perm[i] = b[perm[i]]; });
	Kokkos::fence();

    // Solve permuted system
    solveInPlacePermuted(b_perm);

    // Unpermute solution: x = P^T * x_perm
    Kokkos::parallel_for(
        "b unpermute", Kokkos::RangePolicy<ExecSpace>(0, n), KOKKOS_LAMBDA(const int i) { b[i] = b_perm[perm_inv[i]]; });
	Kokkos::fence();
}

/**
 * Performs forward/backward substitution on permuted system
 * @param b - Permuted right-hand side vector (overwritten with solution)
 */
template <typename T, class MemorySpace>
void SparseLUSolver<T, MemorySpace>::solveInPlacePermuted(const Vector<T, MemorySpace>& b) const
{
    const int n = L_row_ptr_.size() - 1;

    // A loop of size 1 so that calculations are run on GPU
    Kokkos::parallel_for(
        "solveInPlacePermuted", Kokkos::RangePolicy<ExecSpace>(0, 1), KOKKOS_CLASS_LAMBDA(const int) {
            // Forward substitution: L * y = b
            for (int i(0); i < n; ++i) {
                for (int idx = L_row_ptr_[i]; idx < L_row_ptr_[i + 1]; idx++) {
                    b[i] -= L_values_[idx] * b[L_col_idx_[idx]];
                }
            }

            // Backward substitution: U * x = y
            for (int i = n - 1; i >= 0; i--) {
                for (int idx = U_row_ptr_[i]; idx < U_row_ptr_[i + 1]; idx++) {
                    const int col = U_col_idx_[idx];
                    if (col != i) { // Skip diagonal (handled separately)
                        b[i] -= U_values_[idx] * b[col];
                    }
                }
                b[i] /= U_diag_[i]; // Divide by diagonal element
            }
        });
	Kokkos::fence();
}

/**
 * Computes Reverse Cuthill-McKee (RCM) ordering for bandwidth reduction
 * @param A - Input sparse matrix
 * @return Permutation vector
 */
template <typename T, class MemorySpace>
Vector<int, MemorySpace> SparseLUSolver<T, MemorySpace>::computeRCM(const SparseMatrixCSR<T, MemorySpace>& A) const
{
    const int n = A.rows();
    if (n == 0)
        return {};

	auto h_A = A.template mirror_view_and_copy<Kokkos::HostSpace>();

    // Build symmetric adjacency list
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; i++) {
        for (int idx = 0; idx < h_A.row_nz_size(i); idx++) {
            int j = h_A.row_nz_index(i, idx);
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

    Vector<int, Kokkos::HostSpace> h_RCM_order_vec("RCM_order_vec", RCM_order.size());
    for (std::size_t i(0); i < RCM_order.size(); ++i) {
        h_RCM_order_vec[i] = RCM_order[i];
    }

    return Kokkos::create_mirror_view_and_copy(MemorySpace(), h_RCM_order_vec);
}

/**
 * Permutes matrix using RCM ordering (efficient, no sorting)
 * @param A - Original matrix
 * @param perm - Permutation vector
 * @param perm_inv - Inverse permutation
 * @return Permuted CSR matrix
 */
template <typename T, class MemorySpace>
SparseMatrixCSR<T, MemorySpace>
SparseLUSolver<T, MemorySpace>::permuteMatrix(const SparseMatrixCSR<T, MemorySpace>& A,
                                              const Vector<int, MemorySpace>& perm,
                                              const Vector<int, MemorySpace>& perm_inv) const
{
    const int n = A.rows();

    // Compute number of nonzeros per permuted row
    Vector<int, MemorySpace> nz_per_row("nz_per_row", n);
    Kokkos::parallel_for(
        "compute nz_per_row", Kokkos::RangePolicy<ExecSpace>(0, n), KOKKOS_LAMBDA(const int i_new) {
            int i_old         = perm[i_new];
            nz_per_row[i_new] = A.row_nz_size(i_old);
        });
	Kokkos::fence();

    // Construct permuted matrix with preallocated storage
    SparseMatrixCSR<T, MemorySpace> A_perm(n, n, nz_per_row);

    // Fill values and column indices
    Kokkos::parallel_for(
        "SparseLU values and column indices", Kokkos::RangePolicy<ExecSpace>(0, n), KOKKOS_LAMBDA(const int i_new) {
            int i_old = perm[i_new];
            int nnz   = A.row_nz_size(i_old);
            for (int idx = 0; idx < nnz; ++idx) {
                int j_old = A.row_nz_index(i_old, idx);
                T val     = A.row_nz_entry(i_old, idx);
                int j_new = perm_inv[j_old];

                // Find the position in the underlying storage
                A_perm.set_row_nz_entry(i_new, idx, val);
                A_perm.set_row_nz_index(i_new, idx, j_new);
            }
        });
	Kokkos::fence();

    return A_perm;
}

/**
 * Main factorization driver
 * @param A - Permuted matrix to factorize
 */
template <typename T, class MemorySpace>
void SparseLUSolver<T, MemorySpace>::factorize(const SparseMatrixCSR<T, MemorySpace>& A)
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
template <typename T, class MemorySpace>
void SparseLUSolver<T, MemorySpace>::symbolicFactorization(const SparseMatrixCSR<T, MemorySpace>& A,
                                                           std::vector<std::vector<int>>& L_pattern,
                                                           std::vector<std::vector<int>>& U_pattern) const
{
    auto h_A = A.template mirror_view_and_copy<Kokkos::HostSpace>();
    const int n                                  = A.rows();
    L_pattern.resize(n);
    U_pattern.resize(n);

    std::vector<int> marker(n, -1);
    std::vector<int> stk;

    for (int i = 0; i < n; ++i) {
        L_pattern[i].clear();
        U_pattern[i].clear();
        stk.clear();

        // Process original non-zeros in row i
        const int nnz = h_A.row_nz_size(i);
        for (int idx = 0; idx < nnz; ++idx) {
            int col = h_A.row_nz_index(i, idx);
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
template <typename T, class MemorySpace>
void SparseLUSolver<T, MemorySpace>::numericFactorization(const SparseMatrixCSR<T, MemorySpace>& A,
                                                          const std::vector<std::vector<int>>& L_pattern,
                                                          const std::vector<std::vector<int>>& U_pattern)
{
    auto h_A = A.template mirror_view_and_copy<Kokkos::HostSpace>();
    const int n                                  = A.rows();

    // Initialize storage structures
    Vector<int, Kokkos::HostSpace> h_L_row_ptr("L_row_ptr", n + 1);
    Vector<int, Kokkos::HostSpace> h_U_row_ptr("U_row_ptr", n + 1);
    Vector<T, Kokkos::HostSpace> h_U_diag("U_diag", n);

    Kokkos::deep_copy(h_L_row_ptr, 0);
    Kokkos::deep_copy(h_U_row_ptr, 0);
    Kokkos::deep_copy(h_U_diag, 0);

    // Compute row pointers
    for (int i = 0; i < n; i++) {
        h_L_row_ptr[i + 1] = h_L_row_ptr[i] + L_pattern[i].size();
        h_U_row_ptr[i + 1] = h_U_row_ptr[i] + U_pattern[i].size();
    }

    // Allocate memory for values and indices
    Vector<T, Kokkos::HostSpace>   h_L_values  ("L_values", h_L_row_ptr[n]);
    Vector<int, Kokkos::HostSpace> h_L_col_idx ("L_col_idx", h_L_row_ptr[n]);
    Vector<T, Kokkos::HostSpace>   h_U_values  ("U_values", h_U_row_ptr[n]);
    Vector<int, Kokkos::HostSpace> h_U_col_idx ("U_col_idx", h_U_row_ptr[n]);

    // Find start of upper triangular part in U patterns
    Vector<int, Kokkos::HostSpace> U_pattern_start_upper("U_pattern_start_upper", n);
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
        for (int idx = 0; idx < h_A.row_nz_size(i); idx++) {
            int j    = h_A.row_nz_index(i, idx);
            T val    = h_A.row_nz_entry(i, idx);
            dense[j] = val;
            if (marker[j] != i) {
                marker[j] = i;
                indices_used.push(j);
            }
        }

        // Compute L elements
        int L_offset = h_L_row_ptr[i];
        for (const int j : L_pattern[i]) {
            const T Lij = dense[j] / h_U_diag[j];

            h_L_values[L_offset]  = Lij;
            h_L_col_idx[L_offset] = j;
            L_offset++;

            // Update dense row: dense -= Lij * U_row[j] (for columns k > j)
            const int U_update_start_offset = h_U_row_ptr[j] + U_pattern_start_upper[j];
            const int U_row_end_offset      = h_U_row_ptr[j + 1];
            const int update_len            = U_row_end_offset - U_update_start_offset;

            if (update_len <= 0) {
                continue;
            }
            const int* p_U_col = &h_U_col_idx[U_update_start_offset];
            const T* p_U_val   = &h_U_values[U_update_start_offset];

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
        int U_offset        = h_U_row_ptr[i];
        for (int j : U_pattern[i]) {
            T val               = dense[j];
            h_U_values[U_offset]  = val;
            h_U_col_idx[U_offset] = j;
            if (j == i) {
                h_U_diag[i]      = val;
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
        if (std::abs(h_U_diag[i]) < threshold_val) {
            h_U_diag[i] = std::copysign(threshold_val, h_U_diag[i]);
            if (diag_offset != -1) {
                h_U_values[diag_offset] = h_U_diag[i];
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

	L_row_ptr_ = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_L_row_ptr);
	U_row_ptr_ = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_U_row_ptr);
	U_diag_ = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_U_diag);

	L_values_ = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_L_values);
	L_col_idx_ = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_L_col_idx);
	U_values_ = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_U_values);
	U_col_idx_ = Kokkos::create_mirror_view_and_copy(MemorySpace(), h_U_col_idx);

}
} // namespace gmgpolar
