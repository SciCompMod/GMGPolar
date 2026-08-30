#pragma once

/**
 * @brief Batched tridiagonal solver based on a CR-PCR hybrid (Cyclic Reduction
 *        forward/backward sweep around a small Parallel-Cyclic-Reduction core).
 *
 * Motivation (see the accompanying task note / reference paper, "A Hybrid
 * Method for Solving Tridiagonal Systems on the GPU", sections 11.2-11.4.3):
 * pure PCR needs only O(log n) parallel steps but does O(n log n) total work,
 * and its full O(n log n) reduction trajectory (all k1/k2 multipliers, for
 * every step and every equation) is expensive to store. Pure CR is closer to
 * work-efficient (O(n)) but needs 2*log2(n) steps and its later steps under-
 * utilize the vector/SIMD unit because the number of active equations halves
 * every step.
 *
 * CR-PCR runs CR's forward reduction down to a small intermediate system,
 * solves that small system with PCR (its extra work no longer matters once
 * the system is small), then finishes with CR's backward substitution. For
 * this codebase, the decisive benefit is memory: CR's forward reduction only
 * ever updates a shrinking, geometrically-decreasing number of positions per
 * level, and -- summed over all levels -- that is O(n) updates in total, not
 * O(n log n). Storing those updates compactly (one contiguous segment per
 * level, sized to that level's active-position count) instead of in a dense
 * [level][position] grid is what makes this solver's persistent memory
 * footprint O(n), matching BatchedTridiagonalSolverThomas, instead of the
 * O(n log n) footprint of BatchedTridiagonalSolverPCR.
 *
 * ------------------------------------------------------------------------
 * Team sizing: Kokkos::AUTO + strided loops, NOT a fixed team_size (IMPORTANT)
 * ------------------------------------------------------------------------
 * Every phase below (padding/load, CR forward reduction, PCR core, CR
 * backward substitution, write-out) loops over its index range with a
 * strided team loop:
 *
 *     for (int i = team.team_rank(); i < n_padded; i += team.team_size())
 *
 * rather than assuming team_size() == n_padded and using team.team_rank()
 * directly as the position. This is required for portability: host
 * backends (Serial, OpenMP) do not support requesting an arbitrary exact
 * team_size the way CUDA/HIP/SYCL do (an earlier version of this file used
 * TeamPolicy(league_size, n_padded), which aborts at runtime on CPU
 * backends). Kokkos::AUTO lets each backend pick a team_size that fits its
 * execution model -- small (often 1) on CPU, large on GPU -- and the strided
 * loop transparently handles both "one thread owns many positions" (CPU) and
 * "many threads each own one or a few positions" (GPU) without any
 * backend-specific code. This also removes a latent GPU bug: a fixed
 * team_size = n_padded would fail to launch at all once n_padded exceeds the
 * backend's max threads-per-block, which a strided loop has no ceiling on.
 * The PCR core (over the small m-sized intermediate system) uses its own,
 * separate strided loop over t in [0, m) -- it is not tied to the CR phases'
 * loop over i in [0, n_padded), since gather/scatter already bridges the two
 * index spaces explicitly.
 *
 * ------------------------------------------------------------------------
 * Handling of arbitrary matrix_dimension_ (padding, see task section 3):
 * ------------------------------------------------------------------------
 * Classic CR forward-reduction indexing assumes a power-of-two system size.
 * Rather than inventing a bespoke non-power-of-two indexing scheme, every
 * system is padded internally to matrix_dimension_padded_ =
 * next_power_of_two(matrix_dimension_):
 *   - Real equations occupy positions [0, matrix_dimension_).
 *   - Padding positions [matrix_dimension_, matrix_dimension_padded_) are
 *     initialized as *identity* equations: a = 0, b = 1, c = 0 (and d = 0 in
 *     solve()). An identity equation always yields x = 0 and, critically, is
 *     fully decoupled (a = c = 0) from the moment it's created, so it can
 *     never perturb any real equation's reduction no matter how many CR/PCR
 *     steps touch it. This is what makes padding correctness-safe rather
 *     than just convenient.
 *   - The reduction itself runs over all matrix_dimension_padded_ positions
 *     (via the strided team loop described above); only positions with
 *     i < matrix_dimension_ are written back to the caller's rhs.
 *
 * ------------------------------------------------------------------------
 * Compact CR trajectory layout (see task section 4, point 2):
 * ------------------------------------------------------------------------
 * In CR forward reduction, position i is updated in place, and remains
 * "active" (gets re-updated) at level L iff (i+1) is a multiple of
 * 2^(L+1) -- i.e. its 2-adic valuation is >= L+1. A position's *last* update
 * therefore leaves exactly the state backward substitution needs; no
 * per-level snapshot is required, only the last value written to each
 * position (see frozen_sub_lower_/frozen_sub_upper_/frozen_diag_ below).
 * This alone collapses the O(n log n) state storage to O(n).
 *
 * Separately, the *reduction multipliers* (k1, k2) computed at each level ARE
 * needed again later, to replay the same reduction on the right-hand side in
 * solve(). But the number of (level, position) pairs where an update
 * actually happens, summed over all levels, is also O(n) (n/2 + n/4 + ... ~
 * n), because exactly half of the previous level's active positions remain
 * active at the next level. So instead of a dense [level][position] grid
 * (mostly wasted space), the trajectory is stored as one contiguous segment
 * per level, sized to that level's active-position count. cr_level_offsets_
 * holds the prefix-sum start offset of each level's segment (precomputed
 * once in the constructor, identical for every batch entry), so that both
 * setup() and solve() device kernels can independently compute "where does
 * level L's segment start" -- and, using the same i = stride*t + stride-1
 * formula for the t-th active position at level L, always agree on which
 * compact slot corresponds to which system position.
 *
 * ------------------------------------------------------------------------
 * "Frozen" in-place semantics (non-obvious, read before modifying):
 * ------------------------------------------------------------------------
 * After setup() completes, frozen_diag_/frozen_sub_lower_/frozen_sub_upper_
 * do NOT hold "the original matrix" or "the fully-reduced matrix" uniformly
 * -- they hold, *per position*, whatever value that position's forward
 * reduction last computed (or the original loaded value, for positions CR
 * never touched, or the small PCR core's fully-reduced diagonal, for the
 * positions handed off to the PCR sub-solve). Backward substitution in
 * solve() relies on this precisely: see the long comment in solve() for why
 * the *active-position test used for backward substitution is different*
 * from the one used for forward reduction, even though both reuse the same
 * per-level `stride`.
 *
 * ------------------------------------------------------------------------
 * What is unchanged from BatchedTridiagonalSolverPCR (task section 2):
 * ------------------------------------------------------------------------
 *   - The public BatchedTridiagonalSolverBase<T> API and accessor semantics.
 *   - The Sherman-Morrison-Woodbury cyclic wrapper: cyclic systems are still
 *     reduced to two solves of the same plain tridiagonal system (rhs, and
 *     the gamma/corner-derived correction vector), combined via the existing
 *     rank-1 formula. This class only replaces what happens inside "solve a
 *     plain tridiagonal system".
 *   - solve_diagonal() -- pure elementwise scaling, unrelated to any of this.
 *   - The team-per-system parallelization model (Kokkos::TeamPolicy, team
 *     scratch, team_barrier() between read/write phases within a step).
 *
 * @tparam T Scalar type used for matrix coefficients and right-hand sides.
 */

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

namespace gmgpolar
{

/**
 * @brief Rounds n up to the next power of two (returns 1 for n <= 1).
 */
inline int crpcr_next_power_of_two(int n)
{
    if (n <= 1) {
        return 1;
    }
    int p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}

/**
 * @brief Exact base-2 logarithm of a value known to be a power of two
 *        (returns 0 for n <= 1).
 */
inline int crpcr_exact_log2(int n)
{
    int log_value = 0;
    while (n > 1) {
        n >>= 1;
        ++log_value;
    }
    return log_value;
}

/**
 * @brief Clamped left/right neighbor indices, distance `delta` away, used by
 *        the PCR core. Identical convention to BatchedTridiagonalSolverPCR's
 *        pcr_neighbors(): out-of-range neighbors clamp to the array bounds,
 *        which is safe because the corresponding a/c coefficient is zero at
 *        those boundaries.
 */
KOKKOS_INLINE_FUNCTION
void crpcr_pcr_neighbors(int i, int delta, int n, int& iLeft, int& iRight)
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
class BatchedTridiagonalSolverCRPCR : public BatchedTridiagonalSolverBase<T>
{
public:
    /**
     * @param matrix_dimension Real (unpadded) system size, as in the base class.
     * @param batch_count Number of independent systems.
     * @param is_cyclic Whether systems are cyclic (Sherman-Morrison-Woodbury).
     * @param intermediate_system_size Target switch-over size `m` for the CR
     *        forward reduction to hand off to the PCR core (task section 5).
     *        Rounded up to a power of two and clamped to
     *        next_power_of_two(matrix_dimension) if that is smaller -- in
     *        which case CR performs zero levels and this degenerates to pure
     *        PCR on the padded system, which is also the correct behavior
     *        for small matrix_dimension.
     */
    BatchedTridiagonalSolverCRPCR(int matrix_dimension, int batch_count, bool is_cyclic = true,
                                  int intermediate_system_size = 128)
        : BatchedTridiagonalSolverBase<T>(matrix_dimension, batch_count, is_cyclic)
        , matrix_dimension_padded_(crpcr_next_power_of_two(matrix_dimension))
        , intermediate_system_size_(std::min(crpcr_next_power_of_two(std::max(intermediate_system_size, 1)),
                                             crpcr_next_power_of_two(matrix_dimension)))
        , num_cr_levels_(crpcr_exact_log2(crpcr_next_power_of_two(matrix_dimension) /
                                          std::min(crpcr_next_power_of_two(std::max(intermediate_system_size, 1)),
                                                   crpcr_next_power_of_two(matrix_dimension))))
        , num_pcr_steps_(crpcr_exact_log2(std::min(crpcr_next_power_of_two(std::max(intermediate_system_size, 1)),
                                                   crpcr_next_power_of_two(matrix_dimension))))
        , compact_length_(0)
        , cr_level_offsets_("BatchedTridiagonalSolverCRPCR::cr_level_offsets", num_cr_levels_ + 1)
        , frozen_diag_("BatchedTridiagonalSolverCRPCR::frozen_diag",
                       static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(matrix_dimension_padded_))
        , frozen_sub_lower_("BatchedTridiagonalSolverCRPCR::frozen_sub_lower",
                            static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(matrix_dimension_padded_))
        , frozen_sub_upper_("BatchedTridiagonalSolverCRPCR::frozen_sub_upper",
                            static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(matrix_dimension_padded_))
        , gamma_("BatchedTridiagonalSolverCRPCR::gamma", is_cyclic ? batch_count : 0)
        , is_factorized_(false)
    {
        // --- Precompute prefix-sum offsets for the compact CR trajectory. ---
        // offset[0] = 0; offset[L+1] = offset[L] + (active positions at level L).
        // This is host-side, tiny (num_cr_levels_+1 ints), and identical for
        // every batch entry, so both setup() and solve() device kernels can
        // look up "where does level L's compact segment start" via the same
        // small Kokkos::View<int*> without needing an explicit index list.
        std::vector<int> host_offsets(static_cast<std::size_t>(num_cr_levels_) + 1, 0);
        for (int L = 0; L < num_cr_levels_; ++L) {
            const int stride       = 1 << (L + 1);
            const int active_count = matrix_dimension_padded_ / stride;
            host_offsets[L + 1]    = host_offsets[L] + active_count;
        }
        compact_length_ = host_offsets[num_cr_levels_];

        auto offsets_mirror = Kokkos::create_mirror_view(cr_level_offsets_);
        for (int L = 0; L <= num_cr_levels_; ++L) {
            offsets_mirror(L) = host_offsets[L];
        }
        Kokkos::deep_copy(cr_level_offsets_, offsets_mirror);

        // Compact CR multiplier trajectory: O(n) per batch entry (sum over
        // levels of that level's active-position count), not O(n log n).
        cr_k1_trajectory_ =
            Vector<T>("BatchedTridiagonalSolverCRPCR::cr_k1_trajectory",
                      static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(compact_length_));
        cr_k2_trajectory_ =
            Vector<T>("BatchedTridiagonalSolverCRPCR::cr_k2_trajectory",
                      static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(compact_length_));

        // PCR trajectory for the small (size m) intermediate system: O(m log m)
        // per batch entry, negligible next to the O(n) CR trajectory since m
        // is small and fixed (default 32).
        pcr_k1_trajectory_ =
            Vector<T>("BatchedTridiagonalSolverCRPCR::pcr_k1_trajectory",
                      static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(num_pcr_steps_) *
                          static_cast<std::size_t>(intermediate_system_size_));
        pcr_k2_trajectory_ =
            Vector<T>("BatchedTridiagonalSolverCRPCR::pcr_k2_trajectory",
                      static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(num_pcr_steps_) *
                          static_cast<std::size_t>(intermediate_system_size_));

        assign(frozen_diag_, T(0));
        assign(frozen_sub_lower_, T(0));
        assign(frozen_sub_upper_, T(0));
        assign(gamma_, T(0));
        assign(cr_k1_trajectory_, T(0));
        assign(cr_k2_trajectory_, T(0));
        assign(pcr_k1_trajectory_, T(0));
        assign(pcr_k2_trajectory_, T(0));
    }

    /**
     * @brief Performs CR forward reduction down to the intermediate system,
     *        PCR-reduces that intermediate system, and stores everything
     *        needed by solve() to replay the same reduction on a right-hand
     *        side (see the class-level comment for the storage rationale).
     *
     * this->sub_diagonal_ and this->main_diagonal_ (inherited) are read as
     * pure, unmodified input -- exactly as in BatchedTridiagonalSolverPCR,
     * setup() never overwrites them, so setup() may be re-run safely if the
     * matrix changes. (Deviation from the literal task text, which suggested
     * repurposing the inherited main_diagonal_ member to hold the frozen
     * diagonal in place: since main_diagonal_ is populated by callers via
     * BatchedTridiagonalSolverBase's set_main_diagonal(batch_idx, index, ...)
     * accessor, which indexes with the *unpadded* matrix_dimension_ stride,
     * reusing that same storage for the *padded*-length frozen state would
     * require changing accessor semantics, which task section 2 forbids. A
     * dedicated frozen_diag_ member of padded length is used instead; see
     * DEVIATIONS.md.)
     */
    void setup() override
    {
        const int matrix_dimension       = this->matrix_dimension_;
        const int n_padded               = matrix_dimension_padded_;
        const int m                      = intermediate_system_size_;
        const int num_cr_levels          = num_cr_levels_;
        const int num_pcr_steps          = num_pcr_steps_;
        const bool is_cyclic             = this->is_cyclic_;
        const std::size_t compact_length = static_cast<std::size_t>(compact_length_);

        Vector<T> sub_diagonal           = this->sub_diagonal_;
        Vector<T> main_diagonal_in       = this->main_diagonal_;
        Vector<T> frozen_diag            = frozen_diag_;
        Vector<T> frozen_a               = frozen_sub_lower_;
        Vector<T> frozen_c               = frozen_sub_upper_;
        Vector<T> gamma                  = gamma_;
        Vector<T> cr_k1                  = cr_k1_trajectory_;
        Vector<T> cr_k2                  = cr_k2_trajectory_;
        Vector<T> pcr_k1                 = pcr_k1_trajectory_;
        Vector<T> pcr_k2                 = pcr_k2_trajectory_;
        Kokkos::View<int*> level_offsets = cr_level_offsets_;

        // --- Trivial 1x1 systems: nothing to reduce. ---
        if (matrix_dimension == 1) {
            Kokkos::parallel_for(
                "SetupCRPCR_Trivial", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, this->batch_count_),
                KOKKOS_LAMBDA(const int batch_idx) {
                    const int out_offset        = batch_idx * n_padded;
                    frozen_diag(out_offset + 0) = main_diagonal_in(batch_idx * matrix_dimension + 0);
                    frozen_a(out_offset + 0)    = T(0);
                    frozen_c(out_offset + 0)    = T(0);
                    for (int i = 1; i < n_padded; ++i) {
                        // Identity padding, in the unlikely case
                        // matrix_dimension == 1 but intermediate_system_size
                        // rounds n_padded above 1 (it won't, since
                        // next_power_of_two(1) == 1, but keep this robust).
                        frozen_diag(out_offset + i) = T(1);
                        frozen_a(out_offset + i)    = T(0);
                        frozen_c(out_offset + i)    = T(0);
                    }
                });
            Kokkos::fence();
            is_factorized_ = true;
            return;
        }

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;

        // Scratch layout: 6 ping-pong CR arrays (a,b,c x2) of length n_padded,
        // followed by 6 ping-pong PCR-core arrays (a,b,c x2) of length m.
        const std::size_t cr_scratch_elems  = 6ull * static_cast<std::size_t>(n_padded);
        const std::size_t pcr_scratch_elems = 6ull * static_cast<std::size_t>(m);
        const std::size_t scratch_bytes     = (cr_scratch_elems + pcr_scratch_elems) * sizeof(T);

        // Kokkos::AUTO, not a fixed team_size: see the class-level comment on
        // team sizing. Every phase below loops with a strided team loop so it
        // works correctly regardless of what team_size the backend picks.
        TeamPolicy policy(this->batch_count_, Kokkos::AUTO);
        policy.set_scratch_size(1, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

        Kokkos::parallel_for(
            "SetupCRPCR", policy, KOKKOS_LAMBDA(const TeamMember& team) {
                const int batch_idx  = team.league_rank();
                const int rank       = team.team_rank();
                const int team_size  = team.team_size();
                const int sd_offset  = batch_idx * matrix_dimension;
                const int md_offset  = batch_idx * matrix_dimension;
                const int out_offset = batch_idx * n_padded;

                T* scratch = static_cast<T*>(team.team_scratch(1).get_shmem(scratch_bytes));
                T* crA[2]  = {scratch, scratch + 3 * n_padded};
                T* crB[2]  = {scratch + n_padded, scratch + 4 * n_padded};
                T* crC[2]  = {scratch + 2 * n_padded, scratch + 5 * n_padded};
                T* pcrBase = scratch + cr_scratch_elems;
                T* pcrA[2] = {pcrBase, pcrBase + 3 * m};
                T* pcrB[2] = {pcrBase + m, pcrBase + 4 * m};
                T* pcrC[2] = {pcrBase + 2 * m, pcrBase + 5 * m};

                int cur = 0;

                // --- Step 1: load and pad. ---
                // Padding positions are identity equations (a=0,b=1,c=0): they
                // are fully decoupled from creation, so no amount of later
                // reduction can let them perturb a real equation.
                for (int i = rank; i < n_padded; i += team_size) {
                    if (i < matrix_dimension) {
                        crA[cur][i] = (i == 0) ? T(0) : sub_diagonal(sd_offset + i - 1);
                        crC[cur][i] = (i == matrix_dimension - 1) ? T(0) : sub_diagonal(sd_offset + i);
                        crB[cur][i] = main_diagonal_in(md_offset + i);
                    }
                    else {
                        crA[cur][i] = T(0);
                        crB[cur][i] = T(1);
                        crC[cur][i] = T(0);
                    }
                }
                team.team_barrier();

                // --- Step 2: cyclic Sherman-Morrison diagonal adjustment. ---
                // Applied to the *real* boundary positions (0 and
                // matrix_dimension-1), never to the padded system's last
                // position n_padded-1, since padding is not part of the
                // original cyclic system. Done by a single thread (rank==0)
                // since it's two scalar writes into shared scratch.
                if (is_cyclic) {
                    if (rank == 0) {
                        const T cyclic_corner_element = sub_diagonal(sd_offset + matrix_dimension - 1);
                        gamma(batch_idx)              = -main_diagonal_in(md_offset + 0);
                        crB[cur][0] -= gamma(batch_idx);
                        crB[cur][matrix_dimension - 1] -=
                            cyclic_corner_element * cyclic_corner_element / gamma(batch_idx);
                    }
                    team.team_barrier();
                }

                // --- Step 3: CR forward reduction. ---
                // Level L (stride = 2^(L+1)): position i is updated iff
                // (i+1) % stride == 0 (equivalently: 2-adic valuation of
                // (i+1) is >= L+1). A position satisfying this at several
                // consecutive levels gets re-updated each time; the LAST
                // such update is what backward substitution needs, so
                // positions that stop satisfying the condition simply copy
                // their value across unchanged for all remaining levels
                // (see the "else" branch below) -- that's what makes the
                // final buffer contents, after the loop, exactly the frozen
                // state needed (task section 4, point 1).
                for (int L = 0; L < num_cr_levels; ++L) {
                    const int stride = 1 << (L + 1);
                    const int nxt    = 1 - cur;

                    for (int i = rank; i < n_padded; i += team_size) {
                        if ((i + 1) % stride == 0) {
                            const int t = (i + 1) / stride - 1; // compact trajectory slot within this level
                            int iLeft   = i - stride / 2;
                            int iRight  = i + stride / 2;
                            if (iRight > n_padded - 1) {
                                iRight = n_padded - 1;
                            }

                            const T a_i     = crA[cur][i];
                            const T c_i     = crC[cur][i];
                            const T b_left  = crB[cur][iLeft];
                            const T b_right = crB[cur][iRight];

                            const T k1_val = a_i / b_left;
                            const T k2_val = c_i / b_right;

                            const std::size_t traj_idx = static_cast<std::size_t>(batch_idx) * compact_length +
                                                         static_cast<std::size_t>(level_offsets(L)) +
                                                         static_cast<std::size_t>(t);
                            cr_k1(traj_idx)            = k1_val;
                            cr_k2(traj_idx)            = k2_val;

                            crA[nxt][i] = -crA[cur][iLeft] * k1_val;
                            crB[nxt][i] = crB[cur][i] - crC[cur][iLeft] * k1_val - crA[cur][iRight] * k2_val;
                            crC[nxt][i] = -crC[cur][iRight] * k2_val;
                        }
                        else {
                            // Not active at this level: retain whatever value was
                            // last written (do not lose the frozen state).
                            crA[nxt][i] = crA[cur][i];
                            crB[nxt][i] = crB[cur][i];
                            crC[nxt][i] = crC[cur][i];
                        }
                    }

                    team.team_barrier();
                    cur = nxt;
                }

                // --- Step 4: persist frozen forward-reduction state. ---
                // Every position (touched once, many times, or never) now
                // holds its correct frozen a/b/c in the "cur" buffer.
                for (int i = rank; i < n_padded; i += team_size) {
                    frozen_a(out_offset + i)    = crA[cur][i];
                    frozen_diag(out_offset + i) = crB[cur][i];
                    frozen_c(out_offset + i)    = crC[cur][i];
                }
                team.team_barrier();

                // --- Step 5: PCR sub-solve on the size-m intermediate system. ---
                // The m surviving positions are exactly those with valuation
                // >= num_cr_levels, i.e. (i+1) % stride_final == 0 where
                // stride_final = n_padded / m. When num_cr_levels == 0 (CR
                // skipped: n_padded <= requested intermediate size),
                // stride_final == 1 and every position is a "survivor",
                // meaning this step degenerates to plain PCR on the whole
                // padded system, exactly as required.
                //
                // Gather/scatter loop over i in [0, n_padded); the PCR core
                // steps themselves loop separately over t in [0, m) -- a
                // different, independent index space (see class-level note).
                const int stride_final = n_padded / m;
                for (int i = rank; i < n_padded; i += team_size) {
                    if ((i + 1) % stride_final == 0) {
                        const int t = (i + 1) / stride_final - 1;
                        pcrA[0][t]  = crA[cur][i];
                        pcrB[0][t]  = crB[cur][i];
                        pcrC[0][t]  = crC[cur][i];
                    }
                }
                team.team_barrier();

                int pcur = 0;
                for (int step = 0; step < num_pcr_steps; ++step) {
                    const int delta = 1 << step;
                    const int pnxt  = 1 - pcur;

                    for (int t = rank; t < m; t += team_size) {
                        int tLeft, tRight;
                        crpcr_pcr_neighbors(t, delta, m, tLeft, tRight);

                        const T a_cur   = pcrA[pcur][t];
                        const T c_cur   = pcrC[pcur][t];
                        const T b_left  = pcrB[pcur][tLeft];
                        const T b_right = pcrB[pcur][tRight];

                        const T k1_val = a_cur / b_left;
                        const T k2_val = c_cur / b_right;

                        const std::size_t pcr_idx =
                            static_cast<std::size_t>(batch_idx) * static_cast<std::size_t>(num_pcr_steps) *
                                static_cast<std::size_t>(m) +
                            static_cast<std::size_t>(step) * static_cast<std::size_t>(m) + static_cast<std::size_t>(t);
                        pcr_k1(pcr_idx) = k1_val;
                        pcr_k2(pcr_idx) = k2_val;

                        pcrA[pnxt][t] = -pcrA[pcur][tLeft] * k1_val;
                        pcrB[pnxt][t] = pcrB[pcur][t] - pcrC[pcur][tLeft] * k1_val - pcrA[pcur][tRight] * k2_val;
                        pcrC[pnxt][t] = -pcrC[pcur][tRight] * k2_val;
                    }
                    team.team_barrier();
                    pcur = pnxt;
                }

                // Overwrite the CR-frozen diagonal at survivor positions with
                // the fully PCR-reduced diagonal; every other position keeps
                // its CR-only frozen diagonal from step 4.
                for (int i = rank; i < n_padded; i += team_size) {
                    if ((i + 1) % stride_final == 0) {
                        const int t                 = (i + 1) / stride_final - 1;
                        frozen_diag(out_offset + i) = pcrB[pcur][t];
                    }
                }
            });

        Kokkos::fence();
        is_factorized_ = true;
    }

    /**
     * @brief Solves the factored systems for the supplied right-hand side.
     *
     * Mirrors setup()'s CR forward reduction exactly (same per-level active
     * condition, same compact trajectory indices), then replays the PCR
     * core's stored trajectory, then performs CR backward substitution.
     *
     * IMPORTANT: backward substitution's active-position test is NOT the
     * same as forward reduction's, even though both use the same `stride`.
     * See the long comment inside the backward-substitution loop below.
     */
    void solve(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1) override
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        const int effective_batch_count = (this->batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        const int matrix_dimension       = this->matrix_dimension_;
        const int n_padded               = matrix_dimension_padded_;
        const int m                      = intermediate_system_size_;
        const int num_cr_levels          = num_cr_levels_;
        const int num_pcr_steps          = num_pcr_steps_;
        const bool is_cyclic             = this->is_cyclic_;
        const std::size_t compact_length = static_cast<std::size_t>(compact_length_);

        Vector<T> sub_diagonal           = this->sub_diagonal_;
        Vector<T> frozen_diag            = frozen_diag_;
        Vector<T> frozen_a               = frozen_sub_lower_;
        Vector<T> frozen_c               = frozen_sub_upper_;
        Vector<T> gamma                  = gamma_;
        Vector<T> cr_k1                  = cr_k1_trajectory_;
        Vector<T> cr_k2                  = cr_k2_trajectory_;
        Vector<T> pcr_k1                 = pcr_k1_trajectory_;
        Vector<T> pcr_k2                 = pcr_k2_trajectory_;
        Kokkos::View<int*> level_offsets = cr_level_offsets_;

        if (matrix_dimension == 1) {
            Kokkos::parallel_for(
                "SolveCRPCR_Trivial", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, effective_batch_count),
                KOKKOS_LAMBDA(const int k) {
                    const int batch_idx = batch_stride * k + batch_offset;
                    rhs(batch_idx) /= frozen_diag(batch_idx * n_padded + 0);
                });
            Kokkos::fence();
            return;
        }

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;

        if (!is_cyclic) {
            const std::size_t scratch_elems =
                2ull * static_cast<std::size_t>(n_padded) + 2ull * static_cast<std::size_t>(m);
            const std::size_t scratch_bytes = scratch_elems * sizeof(T);

            TeamPolicy policy(effective_batch_count, Kokkos::AUTO);
            policy.set_scratch_size(1, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

            Kokkos::parallel_for(
                "SolveCRPCR_NonCyclic", policy, KOKKOS_LAMBDA(const TeamMember& team) {
                    const int k          = team.league_rank();
                    const int batch_idx  = batch_stride * k + batch_offset;
                    const int rank       = team.team_rank();
                    const int team_size  = team.team_size();
                    const int rhs_offset = batch_idx * matrix_dimension;
                    const int out_offset = batch_idx * n_padded;

                    T* scratch = static_cast<T*>(team.team_scratch(1).get_shmem(scratch_bytes));
                    T* d[2]    = {scratch, scratch + n_padded};
                    T* pd[2]   = {scratch + 2 * n_padded, scratch + 2 * n_padded + m};

                    int cur = 0;
                    for (int i = rank; i < n_padded; i += team_size) {
                        d[cur][i] = (i < matrix_dimension) ? rhs(rhs_offset + i) : T(0);
                    }
                    team.team_barrier();

                    // --- CR forward reduction on d: same levels/order as setup(). ---
                    for (int L = 0; L < num_cr_levels; ++L) {
                        const int stride = 1 << (L + 1);
                        const int nxt    = 1 - cur;

                        for (int i = rank; i < n_padded; i += team_size) {
                            if ((i + 1) % stride == 0) {
                                const int t = (i + 1) / stride - 1;
                                int iLeft   = i - stride / 2;
                                int iRight  = i + stride / 2;
                                if (iRight > n_padded - 1) {
                                    iRight = n_padded - 1;
                                }

                                const std::size_t traj_idx = static_cast<std::size_t>(batch_idx) * compact_length +
                                                             static_cast<std::size_t>(level_offsets(L)) +
                                                             static_cast<std::size_t>(t);
                                const T k1_val             = cr_k1(traj_idx);
                                const T k2_val             = cr_k2(traj_idx);

                                d[nxt][i] = d[cur][i] - d[cur][iLeft] * k1_val - d[cur][iRight] * k2_val;
                            }
                            else {
                                d[nxt][i] = d[cur][i];
                            }
                        }
                        team.team_barrier();
                        cur = nxt;
                    }
                    // d[cur] now holds, for every position, its value as of
                    // its last forward-reduction touch (or the original rhs
                    // value, for positions CR never touched) -- mirroring
                    // exactly what frozen_a/frozen_c/frozen_diag hold.

                    // --- PCR sub-solve on the m surviving positions. ---
                    const int stride_final = n_padded / m;
                    for (int i = rank; i < n_padded; i += team_size) {
                        if ((i + 1) % stride_final == 0) {
                            const int t = (i + 1) / stride_final - 1;
                            pd[0][t]    = d[cur][i];
                        }
                    }
                    team.team_barrier();

                    int pcur = 0;
                    for (int step = 0; step < num_pcr_steps; ++step) {
                        const int delta = 1 << step;
                        const int pnxt  = 1 - pcur;

                        for (int t = rank; t < m; t += team_size) {
                            int tLeft, tRight;
                            crpcr_pcr_neighbors(t, delta, m, tLeft, tRight);

                            const std::size_t pcr_idx = static_cast<std::size_t>(batch_idx) *
                                                            static_cast<std::size_t>(num_pcr_steps) *
                                                            static_cast<std::size_t>(m) +
                                                        static_cast<std::size_t>(step) * static_cast<std::size_t>(m) +
                                                        static_cast<std::size_t>(t);
                            const T k1_val            = pcr_k1(pcr_idx);
                            const T k2_val            = pcr_k2(pcr_idx);

                            pd[pnxt][t] = pd[pcur][t] - pd[pcur][tLeft] * k1_val - pd[pcur][tRight] * k2_val;
                        }
                        team.team_barrier();
                        pcur = pnxt;
                    }

                    // Survivor positions: divide by the fully PCR-reduced
                    // diagonal to get x directly; write it back into the same
                    // working array d[cur] at that position (reused as x from
                    // here on).
                    for (int i = rank; i < n_padded; i += team_size) {
                        if ((i + 1) % stride_final == 0) {
                            const int t = (i + 1) / stride_final - 1;
                            d[cur][i]   = pd[pcur][t] / frozen_diag(out_offset + i);
                        }
                    }
                    team.team_barrier();

                    // --- CR backward substitution. ---
                    //
                    // NOTE (see the class-level derivation): this active-
                    // position test is deliberately DIFFERENT from forward
                    // reduction's, even though the same `stride` is reused.
                    // Forward reduction's test (i+1)%stride==0 selects
                    // positions with valuation(i+1) >= L+1 -- at the deepest
                    // level (L = num_cr_levels-1) this set is *exactly* the
                    // PCR-solved survivor set, which must NOT be recomputed
                    // here. Backward substitution instead solves, at level L,
                    // the positions that were used as (unsolved) NEIGHBORS
                    // during forward level L and then never touched again --
                    // i.e. positions with valuation(i+1) == L exactly, which
                    // is (i+1) % stride == stride/2. Their neighbors at
                    // distance stride/2 have valuation >= L+1, so they are
                    // either PCR survivors or were already solved at a higher
                    // L in this same (decreasing-L) loop -- which is exactly
                    // why the loop must run from num_cr_levels-1 down to 0.
                    for (int L = num_cr_levels - 1; L >= 0; --L) {
                        const int stride = 1 << (L + 1);
                        const int half   = stride / 2;

                        for (int i = rank; i < n_padded; i += team_size) {
                            if ((i + 1) % stride == half) {
                                int pLeft = i - half;
                                if (pLeft < 0) {
                                    pLeft = 0;
                                }
                                int pRight = i + half;
                                if (pRight > n_padded - 1) {
                                    pRight = n_padded - 1;
                                }

                                const T a_i = frozen_a(out_offset + i);
                                const T c_i = frozen_c(out_offset + i);
                                const T b_i = frozen_diag(out_offset + i);

                                d[cur][i] = (d[cur][i] - a_i * d[cur][pLeft] - c_i * d[cur][pRight]) / b_i;
                            }
                        }
                        team.team_barrier();
                    }

                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        rhs(rhs_offset + i) = d[cur][i];
                    }
                });
        }
        else {
            // Cyclic case: run the whole pipeline above twice within one
            // kernel launch, once for the rhs-derived vector and once for the
            // Sherman-Morrison auxiliary vector, sharing the same coefficient
            // trajectories -- exactly as BatchedTridiagonalSolverPCR does.
            const std::size_t scratch_elems =
                4ull * static_cast<std::size_t>(n_padded) + 4ull * static_cast<std::size_t>(m);
            const std::size_t scratch_bytes = scratch_elems * sizeof(T);

            TeamPolicy policy(effective_batch_count, Kokkos::AUTO);
            policy.set_scratch_size(1, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

            Kokkos::parallel_for(
                "SolveCRPCR_Cyclic", policy, KOKKOS_LAMBDA(const TeamMember& team) {
                    const int k          = team.league_rank();
                    const int batch_idx  = batch_stride * k + batch_offset;
                    const int rank       = team.team_rank();
                    const int team_size  = team.team_size();
                    const int rhs_offset = batch_idx * matrix_dimension;
                    const int out_offset = batch_idx * n_padded;

                    const T cyclic_corner_element = sub_diagonal(rhs_offset + matrix_dimension - 1);

                    T* scratch   = static_cast<T*>(team.team_scratch(1).get_shmem(scratch_bytes));
                    T* d_rhs[2]  = {scratch, scratch + 2 * n_padded};
                    T* d_buf[2]  = {scratch + n_padded, scratch + 3 * n_padded};
                    T* pd_rhs[2] = {scratch + 4 * n_padded, scratch + 4 * n_padded + m};
                    T* pd_buf[2] = {scratch + 4 * n_padded + 2 * m, scratch + 4 * n_padded + 3 * m};

                    int cur = 0;
                    for (int i = rank; i < n_padded; i += team_size) {
                        if (i < matrix_dimension) {
                            d_rhs[cur][i] = rhs(rhs_offset + i);
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
                        else {
                            d_rhs[cur][i] = T(0);
                            d_buf[cur][i] = T(0);
                        }
                    }
                    team.team_barrier();

                    // --- CR forward reduction on both vectors. ---
                    for (int L = 0; L < num_cr_levels; ++L) {
                        const int stride = 1 << (L + 1);
                        const int nxt    = 1 - cur;

                        for (int i = rank; i < n_padded; i += team_size) {
                            if ((i + 1) % stride == 0) {
                                const int t = (i + 1) / stride - 1;
                                int iLeft   = i - stride / 2;
                                int iRight  = i + stride / 2;
                                if (iRight > n_padded - 1) {
                                    iRight = n_padded - 1;
                                }

                                const std::size_t traj_idx = static_cast<std::size_t>(batch_idx) * compact_length +
                                                             static_cast<std::size_t>(level_offsets(L)) +
                                                             static_cast<std::size_t>(t);
                                const T k1_val             = cr_k1(traj_idx);
                                const T k2_val             = cr_k2(traj_idx);

                                d_rhs[nxt][i] =
                                    d_rhs[cur][i] - d_rhs[cur][iLeft] * k1_val - d_rhs[cur][iRight] * k2_val;
                                d_buf[nxt][i] =
                                    d_buf[cur][i] - d_buf[cur][iLeft] * k1_val - d_buf[cur][iRight] * k2_val;
                            }
                            else {
                                d_rhs[nxt][i] = d_rhs[cur][i];
                                d_buf[nxt][i] = d_buf[cur][i];
                            }
                        }
                        team.team_barrier();
                        cur = nxt;
                    }

                    // --- PCR sub-solve on both vectors. ---
                    const int stride_final = n_padded / m;
                    for (int i = rank; i < n_padded; i += team_size) {
                        if ((i + 1) % stride_final == 0) {
                            const int t  = (i + 1) / stride_final - 1;
                            pd_rhs[0][t] = d_rhs[cur][i];
                            pd_buf[0][t] = d_buf[cur][i];
                        }
                    }
                    team.team_barrier();

                    int pcur = 0;
                    for (int step = 0; step < num_pcr_steps; ++step) {
                        const int delta = 1 << step;
                        const int pnxt  = 1 - pcur;

                        for (int t = rank; t < m; t += team_size) {
                            int tLeft, tRight;
                            crpcr_pcr_neighbors(t, delta, m, tLeft, tRight);

                            const std::size_t pcr_idx = static_cast<std::size_t>(batch_idx) *
                                                            static_cast<std::size_t>(num_pcr_steps) *
                                                            static_cast<std::size_t>(m) +
                                                        static_cast<std::size_t>(step) * static_cast<std::size_t>(m) +
                                                        static_cast<std::size_t>(t);
                            const T k1_val            = pcr_k1(pcr_idx);
                            const T k2_val            = pcr_k2(pcr_idx);

                            pd_rhs[pnxt][t] =
                                pd_rhs[pcur][t] - pd_rhs[pcur][tLeft] * k1_val - pd_rhs[pcur][tRight] * k2_val;
                            pd_buf[pnxt][t] =
                                pd_buf[pcur][t] - pd_buf[pcur][tLeft] * k1_val - pd_buf[pcur][tRight] * k2_val;
                        }
                        team.team_barrier();
                        pcur = pnxt;
                    }

                    for (int i = rank; i < n_padded; i += team_size) {
                        if ((i + 1) % stride_final == 0) {
                            const int t   = (i + 1) / stride_final - 1;
                            d_rhs[cur][i] = pd_rhs[pcur][t] / frozen_diag(out_offset + i);
                            d_buf[cur][i] = pd_buf[pcur][t] / frozen_diag(out_offset + i);
                        }
                    }
                    team.team_barrier();

                    // --- CR backward substitution on both vectors. --- (see
                    // the detailed note in the non-cyclic path above)
                    for (int L = num_cr_levels - 1; L >= 0; --L) {
                        const int stride = 1 << (L + 1);
                        const int half   = stride / 2;

                        for (int i = rank; i < n_padded; i += team_size) {
                            if ((i + 1) % stride == half) {
                                int pLeft = i - half;
                                if (pLeft < 0) {
                                    pLeft = 0;
                                }
                                int pRight = i + half;
                                if (pRight > n_padded - 1) {
                                    pRight = n_padded - 1;
                                }

                                const T a_i = frozen_a(out_offset + i);
                                const T c_i = frozen_c(out_offset + i);
                                const T b_i = frozen_diag(out_offset + i);

                                d_rhs[cur][i] =
                                    (d_rhs[cur][i] - a_i * d_rhs[cur][pLeft] - c_i * d_rhs[cur][pRight]) / b_i;
                                d_buf[cur][i] =
                                    (d_buf[cur][i] - a_i * d_buf[cur][pLeft] - c_i * d_buf[cur][pRight]) / b_i;
                            }
                        }
                        team.team_barrier();
                    }

                    // --- Sherman-Morrison combination (unchanged logic). ---
                    // Requires x_rhs[0], x_rhs[matrix_dimension-1], x_buf[0],
                    // x_buf[matrix_dimension-1], which are now all available
                    // in d_rhs[cur]/d_buf[cur] regardless of whether those
                    // positions were CR-backward-solved or PCR-core-solved.
                    // Computed redundantly by every thread (cheap scalars,
                    // avoids a broadcast/second barrier).
                    const T dot_product_x_v =
                        d_rhs[cur][0] + cyclic_corner_element / gamma(batch_idx) * d_rhs[cur][matrix_dimension - 1];
                    const T dot_product_u_v =
                        d_buf[cur][0] + cyclic_corner_element / gamma(batch_idx) * d_buf[cur][matrix_dimension - 1];
                    const T factor = dot_product_x_v / (T(1) + dot_product_u_v);

                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        rhs(rhs_offset + i) = d_rhs[cur][i] - factor * d_buf[cur][i];
                    }
                });
        }

        Kokkos::fence();
    }

    /**
     * @brief Solves systems whose matrix has already been reduced to diagonal
     *        form. Unchanged elementwise scaling; only the diagonal storage
     *        it reads from differs (frozen_diag_, at padded stride, instead
     *        of the inherited main_diagonal_).
     */
    void solve_diagonal(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1) override
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        const int effective_batch_count = (this->batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        const int matrix_dimension = this->matrix_dimension_;
        const int n_padded         = matrix_dimension_padded_;
        const bool is_cyclic       = this->is_cyclic_;
        Vector<T> frozen_diag      = frozen_diag_;
        Vector<T> gamma            = gamma_;

        using MDPolicy = Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>;
        MDPolicy policy({0, 0}, {effective_batch_count, matrix_dimension});

        Kokkos::parallel_for(
            "SolveDiagonalCRPCR", policy, KOKKOS_LAMBDA(const int k, const int i) {
                const int batch_idx  = batch_stride * k + batch_offset;
                const int rhs_offset = batch_idx * matrix_dimension;
                const int out_offset = batch_idx * n_padded;

                if (is_cyclic && i == 0) {
                    rhs(rhs_offset) /= frozen_diag(out_offset) + gamma(batch_idx);
                }
                else {
                    rhs(rhs_offset + i) /= frozen_diag(out_offset + i);
                }
            });

        Kokkos::fence();
    }

    // --- Exposed for tests (memory-footprint regression check, etc.) ---
    int matrixDimensionPadded() const
    {
        return matrix_dimension_padded_;
    }
    int intermediateSystemSize() const
    {
        return intermediate_system_size_;
    }
    int numCrLevels() const
    {
        return num_cr_levels_;
    }
    int numPcrSteps() const
    {
        return num_pcr_steps_;
    }
    int compactTrajectoryLength() const
    {
        return compact_length_;
    }

    /**
     * @brief Total bytes of persistent device storage newly introduced by
     *        this class (i.e. not counting the inherited main_diagonal_ /
     *        sub_diagonal_, which are unchanged in size from the base class).
     *        Used by the memory-regression test in task section 11.
     */
    std::size_t persistentDeviceBytes() const
    {
        const std::size_t n_padded  = static_cast<std::size_t>(matrix_dimension_padded_);
        const std::size_t batch     = static_cast<std::size_t>(this->batch_count_);
        const std::size_t m         = static_cast<std::size_t>(intermediate_system_size_);
        const std::size_t pcr_steps = static_cast<std::size_t>(num_pcr_steps_);

        std::size_t bytes = 0;
        bytes += 3 * batch * n_padded * sizeof(T); // frozen_diag_, frozen_sub_lower_, frozen_sub_upper_
        bytes += 2 * batch * static_cast<std::size_t>(compact_length_) * sizeof(T); // cr_k1/k2_trajectory_
        bytes += 2 * batch * pcr_steps * m * sizeof(T); // pcr_k1/k2_trajectory_
        bytes += static_cast<std::size_t>(this->is_cyclic_ ? batch : 0) * sizeof(T); // gamma_
        return bytes;
    }

private:
    int matrix_dimension_padded_; // next_power_of_two(matrix_dimension_)
    int intermediate_system_size_; // m: CR->PCR switch size (power of two, <= matrix_dimension_padded_)
    int num_cr_levels_; // log2(matrix_dimension_padded_ / intermediate_system_size_)
    int num_pcr_steps_; // log2(intermediate_system_size_)
    int compact_length_; // total compact CR trajectory length per batch entry, ~matrix_dimension_padded_

    // Frozen forward-reduction state (see class-level comment): length
    // batch_count_ * matrix_dimension_padded_ each. O(n), not O(n log n).
    Vector<T> frozen_diag_; // "b" at each position, as of its last forward-reduction touch (or PCR-final)
    Vector<T> frozen_sub_lower_; // "a" at each position, as of its last forward-reduction touch
    Vector<T> frozen_sub_upper_; // "c" at each position, as of its last forward-reduction touch

    // Compact CR multiplier trajectory + per-level offsets (see class-level comment). O(n).
    AllocatableVector<T> cr_k1_trajectory_;
    AllocatableVector<T> cr_k2_trajectory_;
    Kokkos::View<int*> cr_level_offsets_; // size num_cr_levels_ + 1, prefix sums, identical across batch entries

    // PCR trajectory for the size-m intermediate system. O(m log m), small and bounded.
    AllocatableVector<T> pcr_k1_trajectory_;
    AllocatableVector<T> pcr_k2_trajectory_;

    Vector<T> gamma_; // Sherman-Morrison diagonal correction, cyclic systems only.
    bool is_factorized_;
};

} // namespace gmgpolar
