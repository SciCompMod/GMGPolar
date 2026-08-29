#pragma once

// -----------------------------------------------------------------------------------------------
// BatchedTridiagonalSolver<T> — Parallel Cyclic Reduction (PCR) implementation.
//
// Drop-in replacement for the previous Thomas-algorithm solver (kept, renamed, in
// batched_tridiagonal_solver_thomas.h as a correctness reference — see the tests in
// test_batched_tridiagonal_solver_pcr.cpp).
//
// Parallelization model
// ----------------------
// The old solver used Kokkos::RangePolicy(0, batch_count_): one thread per system, running the
// fully-sequential Thomas algorithm inside that thread. That is fine when batch_count_ alone gives
// enough independent work to saturate the device, but in this codebase matrix_dimension_ >=
// batch_count_ (few, long lines), so batch_count_ alone leaves the GPU mostly idle.
//
// PCR trades more total work (O(n log n) vs Thomas's O(n)) for parallel depth O(log n) with `n`
// independent workers per step and NO backward-substitution phase. We exploit that by using
// Kokkos::TeamPolicy: one TEAM per system (league_size = batch_count_), and one team THREAD per
// equation within that system (team_size = matrix_dimension_). This recovers the missing
// intra-system parallelism. Do not collapse this back to a RangePolicy-per-system loop that runs
// PCR sequentially inside one thread — that would keep PCR's extra work with none of its extra
// parallelism and would regress versus Thomas.
//
// setup()/solve() split (Strategy A)
// -----------------------------------
// setup() runs the full num_steps_ = ceil(log2(matrix_dimension_)) PCR coefficient reduction on
// a/b/c (never touching the right-hand side) and stores every step's k1/k2 multipliers into
// k1_trajectory_/k2_trajectory_, plus the fully-reduced diagonal into main_diagonal_ (same slot
// that used to hold Thomas's D — same repurposing pattern as before, different contents).
//
// solve() does no coefficient work at all: it replays the stored k1/k2 trajectory to reduce the
// right-hand side d, then divides by the stored diagonal. For the cyclic case, both auxiliary
// solves needed by the Sherman-Morrison reconstruction (rhs and the buffer/gamma vector) are
// driven through the SAME stored trajectory inside a single kernel launch, so the O(n log n)
// coefficient-reduction work happens exactly once per setup() call, never inside solve().
//
// sub_diagonal_ is read-only input from PCR's point of view (unlike the Thomas version, which
// overwrote it in place during factorization) — see the comment on sub_diagonal_ below.
// -----------------------------------------------------------------------------------------------

#include <Kokkos_Core.hpp>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

namespace gmgpolar
{

// ------------------------------------------------------------------------------------------- //
// PCR neighbor helper: clamps i +/- delta into [0, n-1].
//
// This relies on (and inductively preserves, across steps) the convention that the "virtual" a
// at index 0 and c at index n-1 are zero, i.e. equations at the boundary of the system have no
// left/right neighbor respectively. Clamping instead of branching on the boundary avoids
// thread/warp divergence and is the standard published PCR technique; it is also what makes
// non-power-of-two matrix_dimension_ work correctly (see num_steps_ = ceil(log2(n)) below).
// ------------------------------------------------------------------------------------------- //
KOKKOS_INLINE_FUNCTION
void pcr_neighbors(int i, int delta, int n, int& iLeft, int& iRight)
{
    iLeft = i - delta;
    if (iLeft < 0)
        iLeft = 0;
    iRight = i + delta;
    if (iRight >= n)
        iRight = n - 1;
}

template <typename T>
class BatchedTridiagonalSolver
{
public:
    BatchedTridiagonalSolver(int matrix_dimension, int batch_count, bool is_cyclic = true)
        : matrix_dimension_(matrix_dimension)
        , batch_count_(batch_count)
        , main_diagonal_("BatchedTridiagonalSolver::main_diagonal", matrix_dimension * batch_count)
        , sub_diagonal_("BatchedTridiagonalSolver::sub_diagonal", matrix_dimension * batch_count)
        , buffer_("BatchedTridiagonalSolver::buffer", is_cyclic ? matrix_dimension * batch_count : 0)
        , gamma_("BatchedTridiagonalSolver::gamma", is_cyclic ? batch_count : 0)
        , is_cyclic_(is_cyclic)
        , is_factorized_(false)
        , num_steps_(matrix_dimension > 1
                          ? static_cast<int>(std::ceil(std::log2(static_cast<double>(matrix_dimension))))
                          : 0)
        , k1_trajectory_("BatchedTridiagonalSolver::k1_trajectory",
                          static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(num_steps_) *
                              static_cast<std::size_t>(matrix_dimension))
        , k2_trajectory_("BatchedTridiagonalSolver::k2_trajectory",
                          static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(num_steps_) *
                              static_cast<std::size_t>(matrix_dimension))
    {
        assign(main_diagonal_, T(0));
        assign(sub_diagonal_, T(0));
        assign(k1_trajectory_, T(0));
        assign(k2_trajectory_, T(0));

        // ----------------------------------------------------------------------------------- //
        // Team-size limit check (see section 8.1 of the migration plan). PCR's parallelization
        // model requires team_size == matrix_dimension_ (one thread per equation), which may
        // exceed the backend's max threads per team (commonly 1024 on CUDA). We deliberately fail
        // loudly here rather than silently falling back to a strided multi-equation-per-thread
        // scheme, so the limitation is visible instead of silently producing a slow (or, if
        // botched, incorrect) solve. 1024 is a conservative, commonly-safe default; a
        // backend-specific query (e.g. via a throwaway TeamPolicy's team_size_max()) could tighten
        // this further if needed.
        // ----------------------------------------------------------------------------------- //
        constexpr int kMaxTeamSize = 1024;
        if (matrix_dimension_ > kMaxTeamSize) {
            throw std::runtime_error(
                "BatchedTridiagonalSolver: matrix_dimension (" + std::to_string(matrix_dimension_) +
                ") exceeds the maximum supported Kokkos team size (" + std::to_string(kMaxTeamSize) +
                "). PCR's team-per-system parallelization needs one team thread per equation; a "
                "strided multi-equation-per-thread fallback is not implemented.");
        }
    }

    /* ------------------- */
    /* Accessors for sizes */
    /* ------------------- */

    KOKKOS_INLINE_FUNCTION int matrixDimension() const
    {
        return matrix_dimension_;
    }

    KOKKOS_INLINE_FUNCTION int batchCount() const
    {
        return batch_count_;
    }

    /* ---------------------------- */
    /* Accessors for matrix entries */
    /* ---------------------------- */

    KOKKOS_INLINE_FUNCTION const T& main_diagonal(const int batch_idx, const int index) const
    {
        return main_diagonal_(batch_idx * matrix_dimension_ + index);
    }
    KOKKOS_INLINE_FUNCTION void set_main_diagonal(const int batch_idx, const int index, const T& value) const
    {
        main_diagonal_(batch_idx * matrix_dimension_ + index) = value;
    }
    KOKKOS_INLINE_FUNCTION void increase_main_diagonal(const int batch_idx, const int index, const T& value) const
    {
        main_diagonal_(batch_idx * matrix_dimension_ + index) += value;
    }

    // sub_diagonal_ is READ-ONLY input from PCR's point of view. Unlike the old Thomas
    // implementation (which overwrote sub_diagonal_ in place with L's off-diagonal multipliers
    // during setup()), PCR never mutates it: setup() only ever reads a/b/c into team scratch and
    // reduces the scratch copies. This means sub_diagonal_ retains its original, caller-supplied
    // values after setup() — which is required anyway, since re-running setup() after the matrix
    // changes must start from the true matrix, not from a previous factorization's byproduct.
    KOKKOS_INLINE_FUNCTION const T& sub_diagonal(const int batch_idx, const int index) const
    {
        return sub_diagonal_(batch_idx * matrix_dimension_ + index);
    }
    KOKKOS_INLINE_FUNCTION void set_sub_diagonal(const int batch_idx, const int index, const T& value) const
    {
        sub_diagonal_(batch_idx * matrix_dimension_ + index) = value;
    }
    KOKKOS_INLINE_FUNCTION void increase_sub_diagonal(const int batch_idx, const int index, const T& value) const
    {
        sub_diagonal_(batch_idx * matrix_dimension_ + index) += value;
    }

    KOKKOS_INLINE_FUNCTION const T& cyclic_corner(const int batch_idx) const
    {
        return sub_diagonal_(batch_idx * matrix_dimension_ + (matrix_dimension_ - 1));
    }
    KOKKOS_INLINE_FUNCTION T& set_cyclic_corner(const int batch_idx, const T& value) const
    {
        return sub_diagonal_(batch_idx * matrix_dimension_ + (matrix_dimension_ - 1)) = value;
    }
    KOKKOS_INLINE_FUNCTION void increase_cyclic_corner(const int batch_idx, const T& value) const
    {
        sub_diagonal_(batch_idx * matrix_dimension_ + (matrix_dimension_ - 1)) += value;
    }

    // Trajectory accessors, analogous to main_diagonal(batch_idx, index). Mainly useful for
    // tests/debugging; the hot-path kernels below index the flat Views directly (matching the
    // existing code's style of capturing local Vector copies for lambda capture).
    KOKKOS_INLINE_FUNCTION const T& k1(const int batch_idx, const int step, const int index) const
    {
        return k1_trajectory_(static_cast<std::size_t>(batch_idx) * num_steps_ * matrix_dimension_ +
                               static_cast<std::size_t>(step) * matrix_dimension_ + index);
    }
    KOKKOS_INLINE_FUNCTION const T& k2(const int batch_idx, const int step, const int index) const
    {
        return k2_trajectory_(static_cast<std::size_t>(batch_idx) * num_steps_ * matrix_dimension_ +
                               static_cast<std::size_t>(step) * matrix_dimension_ + index);
    }

    /* --------------------------------------------- */
    /* Setup: full PCR coefficient reduction on a/b/c */
    /* --------------------------------------------- */
    // Computes and stores every step's k1/k2 multipliers (Strategy A), and the fully-reduced
    // diagonal into main_diagonal_. For the cyclic case, applies the Sherman-Morrison-Woodbury
    // diagonal adjustment to a *scratch* copy of b before the PCR loop starts — the PCR reduction
    // itself never needs to know it is part of a cyclic system (same separation of concerns as
    // the old Thomas implementation).

    void setup()
    {
        int matrix_dimension = matrix_dimension_;
        int num_steps         = num_steps_;
        bool is_cyclic         = is_cyclic_;

        Vector<T> main_diagonal   = main_diagonal_;
        Vector<T> sub_diagonal    = sub_diagonal_;
        Vector<T> gamma           = gamma_;
        Vector<T> k1_trajectory   = k1_trajectory_;
        Vector<T> k2_trajectory   = k2_trajectory_;

        if (matrix_dimension == 1) {
            // Degenerate case: num_steps_ == 0, nothing to reduce. main_diagonal_ already holds
            // b[0] and solve() will just divide by it. (A cyclic system of dimension 1 is not a
            // meaningful configuration; we intentionally skip the Sherman-Morrison adjustment
            // here rather than apply an ill-defined self-correction.)
            is_factorized_ = true;
            return;
        }

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;

        const int team_size = matrix_dimension;
        // Ping-pong scratch for a, b, c: 2 buffers x 3 arrays x matrix_dimension x sizeof(T).
        const std::size_t scratch_bytes = 2ull * 3ull * static_cast<std::size_t>(matrix_dimension) * sizeof(T);

        TeamPolicy policy(batch_count_, team_size);
        policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

        Kokkos::parallel_for(
            "SetupPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                const int batch_idx = team_member.league_rank();
                const int i         = team_member.team_rank();
                const int offset    = batch_idx * matrix_dimension;

                T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                // Layout: [a0 | b0 | c0 | a1 | b1 | c1], each block matrix_dimension long.
                T* a[2] = {scratch, scratch + 3 * matrix_dimension};
                T* b[2] = {scratch + matrix_dimension, scratch + 4 * matrix_dimension};
                T* c[2] = {scratch + 2 * matrix_dimension, scratch + 5 * matrix_dimension};

                int cur = 0;

                // Load a/b/c from the input matrix. Virtual a[0] and c[n-1] are zero (no
                // left/right neighbor at the boundary), matching the pcr_neighbors() clamping
                // convention.
                a[cur][i] = (i == 0) ? T(0) : sub_diagonal(offset + i - 1);
                b[cur][i] = main_diagonal(offset + i);
                c[cur][i] = (i == matrix_dimension - 1) ? T(0) : sub_diagonal(offset + i);

                team_member.team_barrier();

                if (is_cyclic) {
                    // Sherman-Morrison-Woodbury diagonal adjustment, applied to the *scratch* b
                    // array (not to main_diagonal_ directly) since scratch, not main_diagonal_,
                    // is the live working array during PCR reduction. Same formula as the old
                    // Thomas implementation.
                    if (i == 0) {
                        const T cyclic_corner_element = sub_diagonal(offset + matrix_dimension - 1);
                        gamma(batch_idx)               = -main_diagonal(offset + 0);
                        b[cur][0] -= gamma(batch_idx);
                        b[cur][matrix_dimension - 1] -=
                            cyclic_corner_element * cyclic_corner_element / gamma(batch_idx);
                    }
                    team_member.team_barrier();
                }

                for (int step = 0; step < num_steps; step++) {
                    const int delta = 1 << step;
                    int iLeft, iRight;
                    pcr_neighbors(i, delta, matrix_dimension, iLeft, iRight);

                    const T k1_val = a[cur][i] / b[cur][iLeft];
                    const T k2_val = c[cur][i] / b[cur][iRight];

                    k1_trajectory(static_cast<std::size_t>(batch_idx) * num_steps * matrix_dimension +
                                   static_cast<std::size_t>(step) * matrix_dimension + i) = k1_val;
                    k2_trajectory(static_cast<std::size_t>(batch_idx) * num_steps * matrix_dimension +
                                   static_cast<std::size_t>(step) * matrix_dimension + i) = k2_val;

                    const int nxt = 1 - cur;
                    a[nxt][i]     = -a[cur][iLeft] * k1_val;
                    b[nxt][i]     = b[cur][i] - c[cur][iLeft] * k1_val - a[cur][iRight] * k2_val;
                    c[nxt][i]     = -c[cur][iRight] * k2_val;

                    // PCR reads neighbors' pre-update values within a step, so every thread must
                    // finish writing the "next" buffer before any thread starts reading it as
                    // "current" in the following iteration.
                    team_member.team_barrier();
                    cur = nxt;
                }

                main_diagonal(offset + i) = b[cur][i];
            });

        Kokkos::fence();
        is_factorized_ = true;
    }

    /* -------------------------------------------------------- */
    /* Solve: replay the stored PCR trajectory to reduce the RHS */
    /* -------------------------------------------------------- */
    // No coefficient work here at all — k1/k2 were already computed and stored by setup(). For
    // the cyclic case, the rhs-derived and buffer-derived right-hand sides are pushed through the
    // SAME stored trajectory within one kernel launch, so the O(n log n) reduction only happens
    // once even though two vectors are being solved.

    void solve(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1)
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        const int effective_batch_count = (batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        int matrix_dimension = matrix_dimension_;
        int num_steps         = num_steps_;
        bool is_cyclic         = is_cyclic_;

        Vector<T> main_diagonal = main_diagonal_;
        Vector<T> sub_diagonal  = sub_diagonal_;
        Vector<T> gamma         = gamma_;
        Vector<T> k1_trajectory = k1_trajectory_;
        Vector<T> k2_trajectory = k2_trajectory_;

        if (matrix_dimension == 1) {
            // Degenerate case: num_steps_ == 0, the "solve" is a pure division. Cyclic n=1 is not
            // a meaningful configuration (see setup()); we fall back to the same plain division.
            Kokkos::parallel_for(
                "SolveTrivial", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, effective_batch_count),
                KOKKOS_LAMBDA(const int k) {
                    const int batch_idx = batch_stride * k + batch_offset;
                    rhs(batch_idx) /= main_diagonal(batch_idx);
                });
            Kokkos::fence();
            return;
        }

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;
        const int team_size = matrix_dimension;

        if (!is_cyclic) {
            const std::size_t scratch_bytes = 2ull * static_cast<std::size_t>(matrix_dimension) * sizeof(T);

            TeamPolicy policy(effective_batch_count, team_size);
            policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

            Kokkos::parallel_for(
                "SolveNonCyclicPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                    const int k          = team_member.league_rank();
                    const int batch_idx  = batch_stride * k + batch_offset;
                    const int i          = team_member.team_rank();
                    const int offset     = batch_idx * matrix_dimension;

                    T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                    T* d[2]    = {scratch, scratch + matrix_dimension};

                    int cur   = 0;
                    d[cur][i] = rhs(offset + i);
                    team_member.team_barrier();

                    for (int step = 0; step < num_steps; step++) {
                        const int delta = 1 << step;
                        int iLeft, iRight;
                        pcr_neighbors(i, delta, matrix_dimension, iLeft, iRight);

                        const T k1_val = k1_trajectory(static_cast<std::size_t>(batch_idx) * num_steps *
                                                            matrix_dimension +
                                                        static_cast<std::size_t>(step) * matrix_dimension + i);
                        const T k2_val = k2_trajectory(static_cast<std::size_t>(batch_idx) * num_steps *
                                                            matrix_dimension +
                                                        static_cast<std::size_t>(step) * matrix_dimension + i);

                        const int nxt = 1 - cur;
                        d[nxt][i]     = d[cur][i] - d[cur][iLeft] * k1_val - d[cur][iRight] * k2_val;

                        // PCR reads neighbors' pre-update values within a step; barrier before
                        // swapping which buffer is "current".
                        team_member.team_barrier();
                        cur = nxt;
                    }

                    rhs(offset + i) = d[cur][i] / main_diagonal(offset + i);
                });
        }
        else {
            // Cyclic: reduce both the rhs-derived vector and the buffer/gamma vector through the
            // same stored trajectory in one launch, then perform the Sherman-Morrison
            // reconstruction. cyclic_corner_element is re-read directly from sub_diagonal_ (cheap,
            // and simpler than persisting it as a separate member) since sub_diagonal_ is
            // preserved unmodified by this design.
            const std::size_t scratch_bytes = 4ull * static_cast<std::size_t>(matrix_dimension) * sizeof(T);

            TeamPolicy policy(effective_batch_count, team_size);
            policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

            Kokkos::parallel_for(
                "SolveCyclicPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                    const int k          = team_member.league_rank();
                    const int batch_idx  = batch_stride * k + batch_offset;
                    const int i          = team_member.team_rank();
                    const int offset     = batch_idx * matrix_dimension;

                    T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                    // Layout: [d_rhs(0) | d_buf(0) | d_rhs(1) | d_buf(1)], each block
                    // matrix_dimension long. The "(1)" blocks double as scratch for the final
                    // x_rhs/x_buf values after the reduction loop finishes (see below).
                    T* d_rhs[2] = {scratch, scratch + 2 * matrix_dimension};
                    T* d_buf[2] = {scratch + matrix_dimension, scratch + 3 * matrix_dimension};

                    const T cyclic_corner_element = sub_diagonal(offset + matrix_dimension - 1);

                    int cur = 0;
                    d_rhs[cur][i] = rhs(offset + i);
                    // Initial Sherman-Morrison-Woodbury auxiliary vector: gamma at position 0,
                    // the cyclic corner at position n-1, zero elsewhere. Unlike the sequential
                    // Thomas solve, we do NOT need to hand-chain-multiply this through forward
                    // elimination first — the PCR reduction loop below carries d_buf through the
                    // same steps as d_rhs, which achieves the same effect.
                    if (i == 0)
                        d_buf[cur][i] = gamma(batch_idx);
                    else if (i == matrix_dimension - 1)
                        d_buf[cur][i] = cyclic_corner_element;
                    else
                        d_buf[cur][i] = T(0);

                    team_member.team_barrier();

                    for (int step = 0; step < num_steps; step++) {
                        const int delta = 1 << step;
                        int iLeft, iRight;
                        pcr_neighbors(i, delta, matrix_dimension, iLeft, iRight);

                        const T k1_val = k1_trajectory(static_cast<std::size_t>(batch_idx) * num_steps *
                                                            matrix_dimension +
                                                        static_cast<std::size_t>(step) * matrix_dimension + i);
                        const T k2_val = k2_trajectory(static_cast<std::size_t>(batch_idx) * num_steps *
                                                            matrix_dimension +
                                                        static_cast<std::size_t>(step) * matrix_dimension + i);

                        const int nxt   = 1 - cur;
                        d_rhs[nxt][i] = d_rhs[cur][i] - d_rhs[cur][iLeft] * k1_val - d_rhs[cur][iRight] * k2_val;
                        d_buf[nxt][i] = d_buf[cur][i] - d_buf[cur][iLeft] * k1_val - d_buf[cur][iRight] * k2_val;

                        // Same reasoning as the non-cyclic case: barrier before the buffer swap.
                        team_member.team_barrier();
                        cur = nxt;
                    }

                    const T x_rhs_i = d_rhs[cur][i] / main_diagonal(offset + i);
                    const T x_buf_i = d_buf[cur][i] / main_diagonal(offset + i);

                    // Stash the divided values into the now-unused "other" buffer so every thread
                    // (not just thread 0/n-1) can read x_rhs[0], x_rhs[n-1], x_buf[0], x_buf[n-1]
                    // for the Sherman-Morrison combination below.
                    const int other  = 1 - cur;
                    d_rhs[other][i]  = x_rhs_i;
                    d_buf[other][i]  = x_buf_i;
                    team_member.team_barrier();

                    const T dot_product_x_v = d_rhs[other][0] + cyclic_corner_element / gamma(batch_idx) *
                                                                      d_rhs[other][matrix_dimension - 1];
                    const T dot_product_u_v = d_buf[other][0] + cyclic_corner_element / gamma(batch_idx) *
                                                                      d_buf[other][matrix_dimension - 1];
                    const T factor = dot_product_x_v / (T(1) + dot_product_u_v);

                    rhs(offset + i) = d_rhs[other][i] - factor * d_buf[other][i];
                });
        }
        Kokkos::fence();
    }

    /* ---------------------------- */
    /* Solve: Diagonal Scaling Only */
    /* ---------------------------- */
    // Unchanged from the Thomas version: this path never touches off-diagonal entries or PCR
    // machinery at all, it's pure elementwise scaling, so there's no within-system parallelism to
    // exploit and RangePolicy remains the right tool.

    void solve_diagonal(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1)
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        int effective_batch_count = (batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        int matrix_dimension    = matrix_dimension_;
        Vector<T> main_diagonal = main_diagonal_;
        Vector<T> gamma         = gamma_;

        if (!is_cyclic_) {
            Kokkos::parallel_for(
                "SolveDiagonalNonCyclic", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, effective_batch_count),
                KOKKOS_LAMBDA(const int k) {
                    int batch_idx = batch_stride * k + batch_offset;
                    int offset    = batch_idx * matrix_dimension;
                    for (int i = 0; i < matrix_dimension; i++) {
                        rhs(offset + i) /= main_diagonal(offset + i);
                    }
                });
        }
        else {
            Kokkos::parallel_for(
                "SolveDiagonalCyclic", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, effective_batch_count),
                KOKKOS_LAMBDA(const int k) {
                    int batch_idx = batch_stride * k + batch_offset;
                    int offset    = batch_idx * matrix_dimension;
                    rhs(offset + 0) /= main_diagonal(offset + 0) + gamma(batch_idx);
                    for (int i = 1; i < matrix_dimension; i++) {
                        rhs(offset + i) /= main_diagonal(offset + i);
                    }
                });
        }
        Kokkos::fence();
    }

private:
    int matrix_dimension_;
    int batch_count_;

    Vector<T> main_diagonal_;
    Vector<T> sub_diagonal_;
    Vector<T> buffer_; // kept for API compatibility; unused internally by the PCR solve path,
                        // which builds its buffer/gamma vector directly in team scratch instead.
    Vector<T> gamma_;

    bool is_cyclic_;
    bool is_factorized_;

    int num_steps_; // = ceil(log2(matrix_dimension_)); 0 when matrix_dimension_ == 1.
    Vector<T> k1_trajectory_; // size: batch_count_ * num_steps_ * matrix_dimension_
    Vector<T> k2_trajectory_; // size: batch_count_ * num_steps_ * matrix_dimension_
};
} // namespace gmgpolar
