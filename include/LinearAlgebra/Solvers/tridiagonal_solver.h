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
// Kokkos::TeamPolicy: one TEAM per system (league_size = batch_count_), and team THREADS working
// across the equations within that system. Team size is requested via Kokkos::AUTO rather than
// pinned to matrix_dimension_, since the achievable team size is backend-dependent (large on
// CUDA, often tiny on CPU/OpenMP backends where a "team" maps to hardware threads per core); each
// kernel below distributes the n equations across whatever team size Kokkos actually grants via a
// strided per-thread loop, so a thread simply owns more than one equation when team_size < n. Do
// not collapse this back to a RangePolicy-per-system loop that runs PCR sequentially inside one
// thread — that would keep PCR's extra work with none of its extra parallelism and would regress
// versus Thomas.
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
    if (iLeft < 0) {
        iLeft = 0;
    }

    iRight = i + delta;
    if (iRight >= n) {
        iRight = n - 1;
    }
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
        , gamma_("BatchedTridiagonalSolver::gamma", is_cyclic ? batch_count : 0)
        , is_cyclic_(is_cyclic)
        , is_factorized_(false)
        , num_steps_(
              matrix_dimension > 1 ? static_cast<int>(std::ceil(std::log2(static_cast<double>(matrix_dimension)))) : 0)
        , k1_trajectory_("BatchedTridiagonalSolver::k1_trajectory", static_cast<std::size_t>(batch_count) *
                                                                        static_cast<std::size_t>(num_steps_) *
                                                                        static_cast<std::size_t>(matrix_dimension))
        , k2_trajectory_("BatchedTridiagonalSolver::k2_trajectory", static_cast<std::size_t>(batch_count) *
                                                                        static_cast<std::size_t>(num_steps_) *
                                                                        static_cast<std::size_t>(matrix_dimension))
    {
        assign(main_diagonal_, T(0));
        assign(sub_diagonal_, T(0));
        assign(k1_trajectory_, T(0));
        assign(k2_trajectory_, T(0));

        // ----------------------------------------------------------------------------------- //
        // Team-size handling (section 8.1 of the migration plan).
        //
        // An earlier version of this class launched every TeamPolicy with
        // team_size = matrix_dimension_ and threw at construction time if matrix_dimension_
        // exceeded a hardcoded 1024. That constant was implicitly CUDA-shaped: on CUDA a team
        // does map to a thread block, where team sizes in the hundreds/~1024 are normal. On the
        // OpenMP backend, however, a "team" maps to hardware threads sharing a core, and the
        // actual max team size is typically tiny (single digits) — nothing to do with
        // matrix_dimension_ at all. A fixed constant can therefore be both wrong (too permissive
        // on CPU backends, causing a runtime "Requested too large team size" exception) and
        // pointlessly restrictive (too conservative on GPU backends with a higher real limit).
        //
        // Fix: never assume matrix_dimension_ is an achievable team size. Every TeamPolicy below
        // is launched with Kokkos::AUTO, so Kokkos itself picks a team size the backend can
        // actually satisfy. Inside each kernel, work over the n equations is distributed across
        // whatever team size Kokkos chose via a manual strided loop
        // (`for (i = team_rank(); i < n; i += team_size())`), so a team with only 4 threads (as
        // on your OpenMP run) correctly covers all n equations, each thread just handling
        // ceil(n/team_size) of them sequentially. This works unchanged for n = 1 up to
        // arbitrarily large n, on any backend, with no hardcoded ceiling and no exception.
        // ----------------------------------------------------------------------------------- //
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
        int num_steps        = num_steps_;
        bool is_cyclic       = is_cyclic_;

        Vector<T> main_diagonal = main_diagonal_;
        Vector<T> sub_diagonal  = sub_diagonal_;
        Vector<T> gamma         = gamma_;
        Vector<T> k1_trajectory = k1_trajectory_;
        Vector<T> k2_trajectory = k2_trajectory_;

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

        // Ping-pong scratch for a, b, c: 2 buffers x 3 arrays x matrix_dimension x sizeof(T).
        // Sized by matrix_dimension_ (one slot per equation), independent of however many actual
        // team threads Kokkos ends up giving us below.
        const std::size_t scratch_bytes = 2ull * 3ull * static_cast<std::size_t>(matrix_dimension) * sizeof(T);

        // Kokkos::AUTO lets Kokkos pick a team size the backend can actually satisfy (this can be
        // far smaller than matrix_dimension_, e.g. on CPU/OpenMP backends — see the constructor
        // comment). The strided loops below (`for (i = team_rank(); i < n; i += team_size())`)
        // make the kernel correct for any team size Kokkos chooses, from 1 up to n.
        TeamPolicy policy(batch_count_, Kokkos::AUTO);
        policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

        Kokkos::parallel_for(
            "SetupPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                const int batch_idx = team_member.league_rank();
                const int offset    = batch_idx * matrix_dimension;
                const int team_size = team_member.team_size();
                const int rank      = team_member.team_rank();

                T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                // Layout: [a0 | b0 | c0 | a1 | b1 | c1], each block matrix_dimension long.
                T* a[2] = {scratch, scratch + 3 * matrix_dimension};
                T* b[2] = {scratch + matrix_dimension, scratch + 4 * matrix_dimension};
                T* c[2] = {scratch + 2 * matrix_dimension, scratch + 5 * matrix_dimension};

                int cur = 0;

                // Load a/b/c from the input matrix. Virtual a[0] and c[n-1] are zero (no
                // left/right neighbor at the boundary), matching the pcr_neighbors() clamping
                // convention. Each thread strides over the equations it owns; if team_size < n,
                // a thread simply owns more than one equation and loads them in turn.
                for (int i = rank; i < matrix_dimension; i += team_size) {
                    a[cur][i] = (i == 0) ? T(0) : sub_diagonal(offset + i - 1);
                    b[cur][i] = main_diagonal(offset + i);
                    c[cur][i] = (i == matrix_dimension - 1) ? T(0) : sub_diagonal(offset + i);
                }

                team_member.team_barrier();

                if (is_cyclic) {
                    // Sherman-Morrison-Woodbury diagonal adjustment, applied to the *scratch* b
                    // array (not to main_diagonal_ directly) since scratch, not main_diagonal_,
                    // is the live working array during PCR reduction. Same formula as the old
                    // Thomas implementation. Only equation 0 needs this, so only the thread that
                    // happens to own equation 0 does the work.
                    if (rank == 0) {
                        const T cyclic_corner_element = sub_diagonal(offset + matrix_dimension - 1);
                        gamma(batch_idx)              = -main_diagonal(offset + 0);
                        b[cur][0] -= gamma(batch_idx);
                        b[cur][matrix_dimension - 1] -=
                            cyclic_corner_element * cyclic_corner_element / gamma(batch_idx);
                    }
                    team_member.team_barrier();
                }

                for (int step = 0; step < num_steps; step++) {
                    const int delta = 1 << step;

                    for (int i = rank; i < matrix_dimension; i += team_size) {
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
                    }

                    // PCR reads neighbors' pre-update values within a step, so every thread must
                    // finish writing the "next" buffer (for every equation it owns) before any
                    // thread starts reading it as "current" in the following iteration.
                    team_member.team_barrier();
                    cur = 1 - cur;
                }

                for (int i = rank; i < matrix_dimension; i += team_size) {
                    main_diagonal(offset + i) = b[cur][i];
                }
            });

        Kokkos::fence();
        is_factorized_ = true;
    }

    // void setup()
    // {
    //     int matrix_dimension = matrix_dimension_;
    //     int num_steps        = num_steps_;
    //     bool is_cyclic       = is_cyclic_;

    //     Vector<T> main_diagonal = main_diagonal_;
    //     Vector<T> sub_diagonal  = sub_diagonal_;
    //     Vector<T> gamma         = gamma_;
    //     Vector<T> k1_trajectory = k1_trajectory_;
    //     Vector<T> k2_trajectory = k2_trajectory_;

    //     if (matrix_dimension == 1) {
    //         is_factorized_ = true;
    //         return;
    //     }

    //     using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
    //     using TeamMember = typename TeamPolicy::member_type;

    //     // Symmetric matrix ⇒ PCR reduction preserves a[i] = c[i - delta] (0 if i - delta < 0) at
    //     // every step (proof: substitute the invariant into the newA/newC update formulas — it's
    //     // self-reproducing). So we never need a[] as its own array: it's always a shifted read of
    //     // e[] (= what used to be c[]). Scratch drops from 3 arrays (a,b,c) to 2 (e,b) — 2 buffers
    //     // x 2 arrays x n instead of x 3 — which is a real win since shared memory, not flops, is
    //     // usually what limits team occupancy here. c[]'s own update formula is untouched, so all
    //     // the boundary zero-propagation PCR relies on is unaffected by this change.
    //     const std::size_t scratch_bytes = 2ull * 2ull * static_cast<std::size_t>(matrix_dimension) * sizeof(T);

    //     TeamPolicy policy(batch_count_, Kokkos::AUTO);
    //     policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

    //     Kokkos::parallel_for(
    //         "SetupPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
    //             const int batch_idx = team_member.league_rank();
    //             const int offset    = batch_idx * matrix_dimension;
    //             const int team_size = team_member.team_size();
    //             const int rank      = team_member.team_rank();

    //             T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
    //             // Layout: [e0 | b0 | e1 | b1], each block matrix_dimension long.
    //             T* e[2] = {scratch, scratch + 2 * matrix_dimension};
    //             T* b[2] = {scratch + matrix_dimension, scratch + 3 * matrix_dimension};

    //             int cur = 0;

    //             // e[i] plays the role of the old c[i]; a[i] is never stored — it's read back as
    //             // e[i - delta] inside the step loop below.
    //             for (int i = rank; i < matrix_dimension; i += team_size) {
    //                 e[cur][i] = (i == matrix_dimension - 1) ? T(0) : sub_diagonal(offset + i);
    //                 b[cur][i] = main_diagonal(offset + i);
    //             }

    //             team_member.team_barrier();

    //             if (is_cyclic) {
    //                 if (rank == 0) {
    //                     const T cyclic_corner_element = sub_diagonal(offset + matrix_dimension - 1);
    //                     gamma(batch_idx)              = -main_diagonal(offset + 0);
    //                     b[cur][0] -= gamma(batch_idx);
    //                     b[cur][matrix_dimension - 1] -=
    //                         cyclic_corner_element * cyclic_corner_element / gamma(batch_idx);
    //                 }
    //                 team_member.team_barrier();
    //             }

    //             for (int step = 0; step < num_steps; step++) {
    //                 const int delta = 1 << step;

    //                 for (int i = rank; i < matrix_dimension; i += team_size) {
    //                     int iLeft, iRight;
    //                     pcr_neighbors(i, delta, matrix_dimension, iLeft, iRight);

    //                     // a[i] and a[iRight], derived on the fly from e (symmetry invariant),
    //                     // instead of read from a separately-maintained array.
    //                     const T a_i      = (i >= delta) ? e[cur][i - delta] : T(0);
    //                     const T a_iRight = (iRight >= delta) ? e[cur][iRight - delta] : T(0);
    //                     const T c_i      = e[cur][i];
    //                     const T c_iLeft  = e[cur][iLeft];

    //                     const T k1_val = a_i / b[cur][iLeft];
    //                     const T k2_val = c_i / b[cur][iRight];

    //                     k1_trajectory(static_cast<std::size_t>(batch_idx) * num_steps * matrix_dimension +
    //                                   static_cast<std::size_t>(step) * matrix_dimension + i) = k1_val;
    //                     k2_trajectory(static_cast<std::size_t>(batch_idx) * num_steps * matrix_dimension +
    //                                   static_cast<std::size_t>(step) * matrix_dimension + i) = k2_val;

    //                     const int nxt = 1 - cur;
    //                     // Same formula as the old c[nxt][i] — only a[nxt] has been eliminated.
    //                     e[nxt][i] = -e[cur][iRight] * k2_val;
    //                     b[nxt][i] = b[cur][i] - c_iLeft * k1_val - a_iRight * k2_val;
    //                 }

    //                 team_member.team_barrier();
    //                 cur = 1 - cur;
    //             }

    //             for (int i = rank; i < matrix_dimension; i += team_size) {
    //                 main_diagonal(offset + i) = b[cur][i];
    //             }
    //         });

    //     Kokkos::fence();
    //     is_factorized_ = true;
    // }

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
        int num_steps        = num_steps_;
        bool is_cyclic       = is_cyclic_;

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

        if (!is_cyclic) {
            const std::size_t scratch_bytes = 2ull * static_cast<std::size_t>(matrix_dimension) * sizeof(T);

            // See setup() for why Kokkos::AUTO + strided per-thread loops are used instead of a
            // fixed team_size = matrix_dimension_.
            TeamPolicy policy(effective_batch_count, Kokkos::AUTO);
            policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

            Kokkos::parallel_for(
                "SolveNonCyclicPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                    const int k         = team_member.league_rank();
                    const int batch_idx = batch_stride * k + batch_offset;
                    const int offset    = batch_idx * matrix_dimension;
                    const int team_size = team_member.team_size();
                    const int rank      = team_member.team_rank();

                    T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                    T* d[2]    = {scratch, scratch + matrix_dimension};

                    int cur = 0;
                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        d[cur][i] = rhs(offset + i);
                    }
                    team_member.team_barrier();

                    for (int step = 0; step < num_steps; step++) {
                        const int delta = 1 << step;

                        for (int i = rank; i < matrix_dimension; i += team_size) {
                            int iLeft, iRight;
                            pcr_neighbors(i, delta, matrix_dimension, iLeft, iRight);

                            const T k1_val =
                                k1_trajectory(static_cast<std::size_t>(batch_idx) * num_steps * matrix_dimension +
                                              static_cast<std::size_t>(step) * matrix_dimension + i);
                            const T k2_val =
                                k2_trajectory(static_cast<std::size_t>(batch_idx) * num_steps * matrix_dimension +
                                              static_cast<std::size_t>(step) * matrix_dimension + i);

                            const int nxt = 1 - cur;
                            d[nxt][i]     = d[cur][i] - d[cur][iLeft] * k1_val - d[cur][iRight] * k2_val;
                        }

                        // PCR reads neighbors' pre-update values within a step; barrier before swapping which buffer is "current".
                        team_member.team_barrier();
                        cur = 1 - cur;
                    }

                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        rhs(offset + i) = d[cur][i] / main_diagonal(offset + i);
                    }
                });
        }
        else {
            // Cyclic: reduce both the rhs-derived vector and the buffer/gamma vector through the
            // same stored trajectory in one launch, then perform the Sherman-Morrison
            // reconstruction. cyclic_corner_element is re-read directly from sub_diagonal_ (cheap,
            // and simpler than persisting it as a separate member) since sub_diagonal_ is
            // preserved unmodified by this design.
            const std::size_t scratch_bytes = 4ull * static_cast<std::size_t>(matrix_dimension) * sizeof(T);

            // See setup() for why Kokkos::AUTO + strided per-thread loops are used instead of a
            // fixed team_size = matrix_dimension_.
            TeamPolicy policy(effective_batch_count, Kokkos::AUTO);
            policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

            Kokkos::parallel_for(
                "SolveCyclicPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                    const int k         = team_member.league_rank();
                    const int batch_idx = batch_stride * k + batch_offset;
                    const int offset    = batch_idx * matrix_dimension;
                    const int team_size = team_member.team_size();
                    const int rank      = team_member.team_rank();

                    T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                    // Layout: [d_rhs(0) | d_buf(0) | d_rhs(1) | d_buf(1)], each block
                    // matrix_dimension long. The "(1)" blocks double as scratch for the final
                    // x_rhs/x_buf values after the reduction loop finishes (see below).
                    T* d_rhs[2] = {scratch, scratch + 2 * matrix_dimension};
                    T* d_buf[2] = {scratch + matrix_dimension, scratch + 3 * matrix_dimension};

                    const T cyclic_corner_element = sub_diagonal(offset + matrix_dimension - 1);

                    int cur = 0;
                    // Initial Sherman-Morrison-Woodbury auxiliary vector: gamma at position 0,
                    // the cyclic corner at position n-1, zero elsewhere. Unlike the sequential
                    // Thomas solve, we do NOT need to hand-chain-multiply this through forward
                    // elimination first — the PCR reduction loop below carries d_buf through the
                    // same steps as d_rhs, which achieves the same effect.
                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        d_rhs[cur][i] = rhs(offset + i);
                        if (i == 0) {
                            d_buf[cur][i] = gamma(batch_idx);
                        }
                        else if (i == matrix_dimension - 1) {
                            d_buf[cur][i] = cyclic_corner_element;
                        }
                        else {
                            d_buf[cur][i] = T(0);
                        }
                    }

                    team_member.team_barrier();

                    for (int step = 0; step < num_steps; step++) {
                        const int delta = 1 << step;

                        for (int i = rank; i < matrix_dimension; i += team_size) {
                            int iLeft, iRight;
                            pcr_neighbors(i, delta, matrix_dimension, iLeft, iRight);

                            const T k1_val =
                                k1_trajectory(static_cast<std::size_t>(batch_idx) * num_steps * matrix_dimension +
                                              static_cast<std::size_t>(step) * matrix_dimension + i);
                            const T k2_val =
                                k2_trajectory(static_cast<std::size_t>(batch_idx) * num_steps * matrix_dimension +
                                              static_cast<std::size_t>(step) * matrix_dimension + i);

                            const int nxt = 1 - cur;
                            d_rhs[nxt][i] = d_rhs[cur][i] - d_rhs[cur][iLeft] * k1_val - d_rhs[cur][iRight] * k2_val;
                            d_buf[nxt][i] = d_buf[cur][i] - d_buf[cur][iLeft] * k1_val - d_buf[cur][iRight] * k2_val;
                        }

                        // Same reasoning as the non-cyclic case: barrier before the buffer swap.
                        team_member.team_barrier();
                        cur = 1 - cur;
                    }

                    // Stash the divided values into the now-unused "other" buffer so every thread
                    // (not just whichever thread owns equation 0/n-1) can read x_rhs[0],
                    // x_rhs[n-1], x_buf[0], x_buf[n-1] for the Sherman-Morrison combination below.
                    const int other = 1 - cur;
                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        d_rhs[other][i] = d_rhs[cur][i] / main_diagonal(offset + i);
                        d_buf[other][i] = d_buf[cur][i] / main_diagonal(offset + i);
                    }
                    team_member.team_barrier();

                    const T dot_product_x_v =
                        d_rhs[other][0] + cyclic_corner_element / gamma(batch_idx) * d_rhs[other][matrix_dimension - 1];
                    const T dot_product_u_v =
                        d_buf[other][0] + cyclic_corner_element / gamma(batch_idx) * d_buf[other][matrix_dimension - 1];
                    const T factor = dot_product_x_v / (T(1) + dot_product_u_v);

                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        rhs(offset + i) = d_rhs[other][i] - factor * d_buf[other][i];
                    }
                });
        }
        Kokkos::fence();
    }

    /* ---------------------------- */
    /* Solve: Diagonal Scaling Only */
    /* ---------------------------- */
    // Every (batch, i) scaling is fully independent — no i-1/i+1 coupling at all, unlike the PCR
    // kernels above. That makes the old RangePolicy(0, batch_count_) with a serial inner loop over
    // matrix_dimension_ the worst-parallelized kernel in the file: it exposed only batch_count_
    // independent work items when it could have exposed batch_count_ * matrix_dimension_. Since
    // matrix_dimension_ >= batch_count_ in this codebase (see file header), that serial inner loop
    // was leaving most of the available parallelism on the table. Fixed by flattening to one
    // thread per (batch, i) pair via MDRangePolicy — no team/scratch machinery needed since there's
    // no barrier or neighbor read to synchronize.
    void solve_diagonal(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1)
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        const int effective_batch_count = (batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        const int matrix_dimension = matrix_dimension_;
        const bool is_cyclic       = is_cyclic_;
        Vector<T> main_diagonal    = main_diagonal_;
        Vector<T> gamma            = gamma_;

        using MDPolicy = Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>;
        MDPolicy policy({0, 0}, {effective_batch_count, matrix_dimension});

        Kokkos::parallel_for(
            "SolveDiagonal", policy, KOKKOS_LAMBDA(const int k, const int i) {
                const int batch_idx = batch_stride * k + batch_offset;
                const int offset    = batch_idx * matrix_dimension;

                // Only equation 0 of a cyclic system carries the Sherman-Morrison-Woodbury
                // diagonal correction (gamma); every other (batch, i) pair is a plain divide.
                // This branch is on i itself, so within a given batch's work it diverges for
                // exactly one (batch, 0) thread — negligible, and unavoidable without a second
                // kernel launch just to special-case a single index per system.
                if (is_cyclic && i == 0) {
                    rhs(offset) /= main_diagonal(offset) + gamma(batch_idx);
                }
                else {
                    rhs(offset + i) /= main_diagonal(offset + i);
                }
            });

        Kokkos::fence();
    }

private:
    int matrix_dimension_;
    int batch_count_;

    Vector<T> main_diagonal_;
    Vector<T> sub_diagonal_;

    Vector<T> gamma_;

    bool is_cyclic_;
    bool is_factorized_;

    int num_steps_; // = ceil(log2(matrix_dimension_)); 0 when matrix_dimension_ == 1.
    Vector<T> k1_trajectory_; // size: batch_count_ * num_steps_ * matrix_dimension_
    Vector<T> k2_trajectory_; // size: batch_count_ * num_steps_ * matrix_dimension_
};

} // namespace gmgpolar
