#pragma once

/**
 * @brief Batched tridiagonal solver based on a Cyclic-Reduction / Parallel-Cyclic-Reduction
 *        (CR-PCR) hybrid.
 *
 * This solver replaces the O(n log n)-per-system storage of full Parallel Cyclic Reduction
 * (see BatchedTridiagonalSolverPCR) with an O(n)-per-system footprint, while keeping most of
 * PCR's parallel-depth advantage over sequential Thomas elimination.
 *
 * The approach follows "A Hybrid Method for Solving Tridiagonal Systems on the GPU"
 * (sections 11.2, 11.3.4, 11.4.3): classic Cyclic Reduction (CR) forward-reduces a system of
 * size n down to a small intermediate system of size m (a tunable constructor parameter,
 * `intermediate_system_size`), which is then solved directly with PCR (whose O(m log m) work
 * is cheap because m is small), followed by CR backward substitution to recover the full
 * solution.
 *
 * Why this collapses storage from O(n log n) to O(n): unlike PCR, where *every* equation is
 * updated at *every* one of the O(log n) steps (hence O(n log n) stored multipliers), CR's
 * forward reduction only ever updates a shrinking, geometrically-decreasing number of
 * equations per level (n/2 at level 0, n/4 at level 1, ...). Summed over all levels this is
 * O(n), not O(n log n) -- provided the multiplier trajectory is stored *compactly* (one
 * contiguous segment per level, sized to that level's active-equation count) rather than as a
 * dense [level][equation] grid. See the "Compact CR trajectory storage" comment below for the
 * exact layout.
 *
 * ---------------------------------------------------------------------------------------------
 * IMPORTANT: this implementation deliberately corrects one detail relative to a literal reading
 * of the original task write-up, and the correction is significant enough to explain here as
 * well as at the point of use (see "CR backward substitution" in solve()):
 *
 * The forward-reduction step at level L eliminates positions i-stride/2 and i+stride/2 (the
 * active positions of the *previous*, finer level -- or the original equations, for L == 0)
 * from position i's equation. Algebraically, this makes position i's *new* equation relate
 * x[i-stride] and x[i+stride] -- the FULL current-level stride, not stride/2 -- because
 * eliminating position (i - stride/2)'s equation pulls in *its* neighbor at
 * (i - stride/2) - stride/2 = i - stride (and symmetrically on the right). This is visible
 * directly in the update itself: `a[next][i] = -a[cur][iLeft] * k1` is literally the
 * coefficient of x at position (iLeft - stride/2) = (i - stride), inherited from
 * position iLeft's own equation.
 *
 * Consequently, CR backward substitution must reconstruct x[i] using neighbors at
 * i -+ stride (the full stride of the level where i was *last* updated), not i -+ stride/2.
 * Reusing the forward step's i -+ stride/2 pair for backward substitution -- as a literal
 * reading of "same iLeft/iRight as forward reduction level L" might suggest -- would relate
 * x[i] to neighbors that generally have *not* been resolved yet at that point in the backward
 * sweep (they belong to a finer, not-yet-processed level), producing incorrect results.
 *
 * A second, related correction: a position i that satisfies level L's "active position" index
 * formula is not necessarily *last touched* at level L -- it may be touched again (overwritten)
 * at a coarser level L' > L, all the way up to the CR/PCR handoff level. Backward substitution
 * must therefore visit each position exactly once, at its true last-touch level, which this
 * implementation identifies via a simple parity test (see solve() for the derivation and
 * comments). Positions that survive all the way to the CR/PCR handoff (the `m` positions PCR
 * solves directly) require no CR backward step at all -- PCR's full reduction already gives
 * their exact solution.
 *
 * This was verified by hand against a concrete n = 8 example before being implemented; see the
 * accompanying deviations note for the full derivation.
 * ---------------------------------------------------------------------------------------------
 *
 * What is unchanged from BatchedTridiagonalSolverPCR:
 *  - The team-per-system / thread-per-equation Kokkos::TeamPolicy parallelization model, team
 *    scratch memory, and team_barrier()-separated read/write phases within a step.
 *  - The Sherman-Morrison-Woodbury handling of cyclic systems (reducing a cyclic system to two
 *    solves of the same plain tridiagonal system, combined via the existing rank-1 formula).
 *    CR-PCR only changes what happens inside "solve a plain tridiagonal system".
 *  - solve_diagonal(): pure elementwise scaling, unrelated to the reduction machinery.
 *
 * @tparam T Scalar type used for matrix coefficients and right-hand sides.
 */

#include <Kokkos_Core.hpp>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

// Reused for pcr_neighbors() in the small intermediate-system PCR sub-solve, and because
// BatchedTridiagonalSolverPCR is kept around as a correctness/performance reference
// (see task deliverable notes) until this implementation is validated in production.
#include "tridiagonal_solver_pcr.h"

namespace gmgpolar
{

/**
 * @brief Rounds n up to the next power of two (returns 1 for n <= 1).
 *
 * Host-only helper used once per solver instance, in the constructor, to determine the padded
 * system size (see the class-level comment and section 3 of the design note for why padding to
 * a power of two is necessary and safe).
 */
inline int next_power_of_two(int n)
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
 * @brief Returns whether position i is one of the "active" (updated) positions at CR level L.
 *
 * Level L (0-indexed, stride = 2^(L+1)) has n_padded / stride active positions, located at
 * i = stride*t + stride - 1 for t = 0 .. active_count-1. This tests the inverse: given i,
 * whether it is expressible in that form for the supplied stride/active_count.
 */
KOKKOS_INLINE_FUNCTION
bool cr_is_active_at_level(int i, int stride, int active_count)
{
    if (((i + 1) % stride) != 0) {
        return false;
    }
    const int t = (i + 1) / stride - 1;
    return (t >= 0 && t < active_count);
}

template <typename T>
class BatchedTridiagonalSolverCRPCR : public BatchedTridiagonalSolverBase<T>
{
public:
    /**
     * @param matrix_dimension Size of each (unpadded) tridiagonal system.
     * @param batch_count Number of independent systems.
     * @param is_cyclic Whether the systems are cyclic (Sherman-Morrison-Woodbury handling).
     * @param intermediate_system_size The size `m` of the intermediate system handed off from
     *        CR to PCR. Must be a power of two. Default 32 (see design note section 5 for the
     *        reasoning behind this default and how to tune it). If the padded matrix dimension
     *        is already <= intermediate_system_size, CR is skipped entirely and the solver
     *        degrades gracefully to pure PCR on the padded system.
     */
    BatchedTridiagonalSolverCRPCR(int matrix_dimension, int batch_count, bool is_cyclic = true,
                                  int intermediate_system_size = 32)
        : BatchedTridiagonalSolverBase<T>(matrix_dimension, batch_count, is_cyclic)
        , matrix_dimension_padded_(next_power_of_two(matrix_dimension))
        , intermediate_system_size_(intermediate_system_size)
        , num_cr_levels_(0)
        , num_pcr_steps_(0)
        , cr_trajectory_length_per_batch_(0)
        , gamma_("BatchedTridiagonalSolverCRPCR::gamma", is_cyclic ? batch_count : 0)
        , is_factorized_(false)
    {
        if (intermediate_system_size_ <= 0 || (intermediate_system_size_ & (intermediate_system_size_ - 1)) != 0) {
            throw std::invalid_argument("intermediate_system_size must be a positive power of two");
        }

        // If the padded system is already no larger than the requested intermediate size,
        // skip CR entirely: run pure PCR directly on the padded system. This also correctly
        // and uniformly handles matrix_dimension == 1 and 2 without special-casing elsewhere.
        if (matrix_dimension_padded_ <= intermediate_system_size_) {
            intermediate_system_size_ = matrix_dimension_padded_;
            num_cr_levels_            = 0;
        }
        else {
            num_cr_levels_ = static_cast<int>(
                std::round(std::log2(static_cast<double>(matrix_dimension_padded_) / intermediate_system_size_)));
        }

        num_pcr_steps_ = (intermediate_system_size_ > 1)
                             ? static_cast<int>(std::round(std::log2(static_cast<double>(intermediate_system_size_))))
                             : 0;

        // Precompute prefix-sum offsets of the per-level active-position counts, so that
        // setup() and solve() can each independently compute "where does level L's segment
        // start in the compact trajectory array" -- see the class-level comment on why this
        // compact (not dense [level][position]) layout is what makes storage O(n).
        std::vector<int> offsets_host(num_cr_levels_ + 1, 0);
        for (int L = 0; L < num_cr_levels_; ++L) {
            const int stride       = 1 << (L + 1);
            const int active_count = matrix_dimension_padded_ / stride;
            offsets_host[L + 1]    = offsets_host[L] + active_count;
        }
        cr_trajectory_length_per_batch_ = offsets_host[num_cr_levels_]; // ~ matrix_dimension_padded_

        cr_level_offsets_ = Kokkos::View<int*>("BatchedTridiagonalSolverCRPCR::cr_level_offsets", num_cr_levels_ + 1);
        auto host_offsets = Kokkos::create_mirror_view(cr_level_offsets_);
        for (int L = 0; L <= num_cr_levels_; ++L) {
            host_offsets(L) = offsets_host[L];
        }
        Kokkos::deep_copy(cr_level_offsets_, host_offsets);

        const std::size_t padded_len =
            static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(matrix_dimension_padded_);
        const std::size_t cr_traj_len =
            static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(cr_trajectory_length_per_batch_);
        const std::size_t pcr_traj_len = static_cast<std::size_t>(batch_count) *
                                         static_cast<std::size_t>(num_pcr_steps_) *
                                         static_cast<std::size_t>(intermediate_system_size_);

        // Frozen forward-reduction state (O(n)): the last value written to each padded
        // position's a/b/c during CR forward reduction (or the initial, untouched value, for
        // positions CR never updated). This is all backward substitution needs -- no
        // per-level snapshot required.
        frozen_sub_lower_     = Vector<T>("BatchedTridiagonalSolverCRPCR::frozen_sub_lower", padded_len);
        frozen_main_diagonal_ = Vector<T>("BatchedTridiagonalSolverCRPCR::frozen_main_diagonal", padded_len);
        frozen_sub_upper_     = Vector<T>("BatchedTridiagonalSolverCRPCR::frozen_sub_upper", padded_len);
        assign(frozen_sub_lower_, T(0));
        assign(frozen_main_diagonal_, T(0));
        assign(frozen_sub_upper_, T(0));

        // Compact CR multiplier trajectory (O(n), see class-level comment).
        cr_k1_trajectory_ = Vector<T>("BatchedTridiagonalSolverCRPCR::cr_k1_trajectory", cr_traj_len);
        cr_k2_trajectory_ = Vector<T>("BatchedTridiagonalSolverCRPCR::cr_k2_trajectory", cr_traj_len);
        assign(cr_k1_trajectory_, T(0));
        assign(cr_k2_trajectory_, T(0));

        // PCR trajectory for the small (size m) intermediate system: O(m log m), negligible
        // next to the O(n) CR trajectory since m is small and fixed.
        pcr_k1_trajectory_ = Vector<T>("BatchedTridiagonalSolverCRPCR::pcr_k1_trajectory", pcr_traj_len);
        pcr_k2_trajectory_ = Vector<T>("BatchedTridiagonalSolverCRPCR::pcr_k2_trajectory", pcr_traj_len);
        assign(pcr_k1_trajectory_, T(0));
        assign(pcr_k2_trajectory_, T(0));

        assign(gamma_, T(0));
    }

    /**
     * @brief Performs CR forward reduction down to the intermediate system, then PCR-reduces
     *        that intermediate system, storing everything needed to replay both stages in
     *        solve().
     *
     * NOTE: unlike BatchedTridiagonalSolverPCR, this implementation does *not* repurpose the
     * inherited main_diagonal_/sub_diagonal_ members to hold padded/frozen state. Doing so
     * would change their effective size and per-batch stride from what
     * BatchedTridiagonalSolverBase's public accessors (main_diagonal(batch_idx, index), etc.)
     * assume, which are used by external code to populate the matrix before setup() is called --
     * silently breaking that population step. Instead, main_diagonal_/sub_diagonal_ are used
     * exactly as the base class documents (pure input, read here, never resized or overwritten
     * with padded/reduced content), and a dedicated set of padded-size members
     * (frozen_main_diagonal_, frozen_sub_lower_, frozen_sub_upper_) holds the frozen state. See
     * the deviations note for more detail; this keeps the base class's accessor semantics
     * exactly as specified (a non-negotiable requirement), at the cost of one additional O(n)
     * array (frozen_main_diagonal_) that a literal member-reuse would have avoided -- a small
     * price for not silently corrupting external population of the matrix.
     */
    void setup() override
    {
        const int matrix_dimension = this->matrix_dimension_;
        const int n_padded         = matrix_dimension_padded_;
        const int m                = intermediate_system_size_;
        const int num_cr_levels    = num_cr_levels_;
        const int num_pcr_steps    = num_pcr_steps_;
        const int cr_traj_len      = cr_trajectory_length_per_batch_;
        const bool is_cyclic       = this->is_cyclic_;

        Vector<T> main_diagonal_in       = this->main_diagonal_;
        Vector<T> sub_diagonal_in        = this->sub_diagonal_;
        Vector<T> frozen_a               = frozen_sub_lower_;
        Vector<T> frozen_b               = frozen_main_diagonal_;
        Vector<T> frozen_c               = frozen_sub_upper_;
        Vector<T> cr_k1                  = cr_k1_trajectory_;
        Vector<T> cr_k2                  = cr_k2_trajectory_;
        Vector<T> pcr_k1                 = pcr_k1_trajectory_;
        Vector<T> pcr_k2                 = pcr_k2_trajectory_;
        Kokkos::View<int*> level_offsets = cr_level_offsets_;
        Vector<T> gamma                  = gamma_;

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;

        // Team scratch: ping-pong a,b,c arrays of length n_padded for CR forward reduction,
        // plus a second, smaller ping-pong a,b,c region of length m for the PCR sub-solve on
        // the intermediate system.
        const std::size_t main_scratch_bytes = 6ull * static_cast<std::size_t>(n_padded) * sizeof(T);
        const std::size_t pcr_scratch_bytes  = 6ull * static_cast<std::size_t>(m) * sizeof(T);
        const std::size_t scratch_bytes      = main_scratch_bytes + pcr_scratch_bytes;

        // Team size is the padded dimension (section 3): every padded position, real or
        // identity, participates in the reduction, avoiding divergent indexing logic. Loops
        // below are still strided by the *actual* team size in case a backend caps it below
        // the requested value.
        TeamPolicy policy(this->batch_count_, n_padded);
        policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

        Kokkos::parallel_for(
            "SetupCRPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                const int batch_idx     = team_member.league_rank();
                const int rank          = team_member.team_rank();
                const int team_size     = team_member.team_size();
                const int in_offset     = batch_idx * matrix_dimension;
                const int padded_offset = batch_idx * n_padded;

                T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                // Main ping-pong buffers, each length n_padded: [a0|b0|c0 | a1|b1|c1].
                T* a[2] = {scratch, scratch + 3 * n_padded};
                T* b[2] = {scratch + n_padded, scratch + 4 * n_padded};
                T* c[2] = {scratch + 2 * n_padded, scratch + 5 * n_padded};
                // Intermediate-system (size m) ping-pong buffers, placed right after.
                T* pcr_base = scratch + 6 * n_padded;
                T* pa[2]    = {pcr_base, pcr_base + 3 * m};
                T* pb[2]    = {pcr_base + m, pcr_base + 4 * m};
                T* pc[2]    = {pcr_base + 2 * m, pcr_base + 5 * m};

                int cur = 0;

                // --- Step 1: load and pad. ---
                // Real equations occupy [0, matrix_dimension). Padding positions
                // [matrix_dimension, n_padded) are identity equations (a=0, b=1, c=0): fully
                // decoupled from the moment they're created, so no number of CR/PCR steps can
                // let them perturb any real equation's reduction. This is what makes padding
                // to a power of two safe.
                for (int i = rank; i < n_padded; i += team_size) {
                    if (i < matrix_dimension) {
                        a[cur][i] = (i == 0) ? T(0) : sub_diagonal_in(in_offset + i - 1);
                        b[cur][i] = main_diagonal_in(in_offset + i);
                        c[cur][i] = (i == matrix_dimension - 1) ? T(0) : sub_diagonal_in(in_offset + i);
                    }
                    else {
                        a[cur][i] = T(0);
                        b[cur][i] = T(1);
                        c[cur][i] = T(0);
                    }
                }
                team_member.team_barrier();

                // --- Step 2: cyclic Sherman-Morrison diagonal adjustment. ---
                // Applied to the *real* boundary positions (matrix_dimension - 1), not the
                // padded system's last position -- padding positions are not part of the
                // original cyclic system.
                if (is_cyclic) {
                    if (rank == 0) {
                        const T cyclic_corner_element = sub_diagonal_in(in_offset + matrix_dimension - 1);
                        gamma(batch_idx)              = -main_diagonal_in(in_offset + 0);
                        b[cur][0] -= gamma(batch_idx);
                        b[cur][matrix_dimension - 1] -=
                            cyclic_corner_element * cyclic_corner_element / gamma(batch_idx);
                    }
                    team_member.team_barrier();
                }

                // --- Step 3: CR forward reduction. ---
                for (int L = 0; L < num_cr_levels; ++L) {
                    const int stride       = 1 << (L + 1);
                    const int active_count = n_padded / stride;
                    const int nxt          = 1 - cur;

                    for (int t = rank; t < active_count; t += team_size) {
                        const int i     = stride * t + stride - 1;
                        const int iLeft = i - stride / 2;
                        int iRight      = i + stride / 2;
                        if (iRight > n_padded - 1) {
                            iRight = n_padded - 1;
                        }

                        const T k1_val = a[cur][i] / b[cur][iLeft];
                        const T k2_val = c[cur][i] / b[cur][iRight];

                        cr_k1(static_cast<std::size_t>(batch_idx) * cr_traj_len + level_offsets(L) + t) = k1_val;
                        cr_k2(static_cast<std::size_t>(batch_idx) * cr_traj_len + level_offsets(L) + t) = k2_val;

                        a[nxt][i] = -a[cur][iLeft] * k1_val;
                        b[nxt][i] = b[cur][i] - c[cur][iLeft] * k1_val - a[cur][iRight] * k2_val;
                        c[nxt][i] = -c[cur][iRight] * k2_val;
                    }

                    // Positions not active at this level keep whatever value they last had.
                    for (int i = rank; i < n_padded; i += team_size) {
                        if (!cr_is_active_at_level(i, stride, active_count)) {
                            a[nxt][i] = a[cur][i];
                            b[nxt][i] = b[cur][i];
                            c[nxt][i] = c[cur][i];
                        }
                    }

                    team_member.team_barrier();
                    cur = nxt;
                }

                // --- Step 4: freeze forward-reduction state. ---
                // scratch now holds, for every padded position, its value as of its last
                // forward-reduction touch (or its initial load, if CR never touched it).
                for (int i = rank; i < n_padded; i += team_size) {
                    frozen_a(padded_offset + i) = a[cur][i];
                    frozen_b(padded_offset + i) = b[cur][i];
                    frozen_c(padded_offset + i) = c[cur][i];
                }
                team_member.team_barrier();

                // --- Step 5: PCR sub-solve setup on the size-m intermediate system. ---
                const int stride_final = n_padded / m;

                for (int t = rank; t < m; t += team_size) {
                    const int i = stride_final * t + stride_final - 1;
                    pa[0][t]    = a[cur][i];
                    pb[0][t]    = b[cur][i];
                    pc[0][t]    = c[cur][i];
                }
                team_member.team_barrier();

                int pcur = 0;
                for (int step = 0; step < num_pcr_steps; ++step) {
                    const int delta = 1 << step;
                    const int pnxt  = 1 - pcur;

                    for (int t = rank; t < m; t += team_size) {
                        int tLeft, tRight;
                        pcr_neighbors(t, delta, m, tLeft, tRight);

                        const T a_t      = (t >= delta) ? pa[pcur][t - delta] : T(0);
                        const T a_tRight = (tRight >= delta) ? pa[pcur][tRight - delta] : T(0);
                        const T c_t      = pc[pcur][t];
                        const T c_tLeft  = pc[pcur][tLeft];

                        const T k1_val = a_t / pb[pcur][tLeft];
                        const T k2_val = c_t / pb[pcur][tRight];

                        pcr_k1(static_cast<std::size_t>(batch_idx) * num_pcr_steps * m +
                               static_cast<std::size_t>(step) * m + t) = k1_val;
                        pcr_k2(static_cast<std::size_t>(batch_idx) * num_pcr_steps * m +
                               static_cast<std::size_t>(step) * m + t) = k2_val;

                        pa[pnxt][t] = -pa[pcur][tLeft] * k1_val;
                        pb[pnxt][t] = pb[pcur][t] - c_tLeft * k1_val - a_tRight * k2_val;
                        pc[pnxt][t] = -pc[pcur][tRight] * k2_val;
                    }

                    team_member.team_barrier();
                    pcur = pnxt;
                }

                // Fully-PCR-reduced diagonal at the m surviving positions overwrites the
                // CR-frozen b there: those positions are now completely reduced.
                for (int t = rank; t < m; t += team_size) {
                    const int i                 = stride_final * t + stride_final - 1;
                    frozen_b(padded_offset + i) = pb[pcur][t];
                }
            });

        Kokkos::fence();
        is_factorized_ = true;
    }

    /**
     * @brief Solves the factored systems for the supplied right-hand side.
     */
    void solve(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1) override
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        const int effective_batch_count = (this->batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        const int matrix_dimension = this->matrix_dimension_;
        const int n_padded         = matrix_dimension_padded_;
        const int m                = intermediate_system_size_;
        const int num_cr_levels    = num_cr_levels_;
        const int num_pcr_steps    = num_pcr_steps_;
        const int cr_traj_len      = cr_trajectory_length_per_batch_;
        const bool is_cyclic       = this->is_cyclic_;

        Vector<T> frozen_a               = frozen_sub_lower_;
        Vector<T> frozen_b               = frozen_main_diagonal_;
        Vector<T> frozen_c               = frozen_sub_upper_;
        Vector<T> cr_k1                  = cr_k1_trajectory_;
        Vector<T> cr_k2                  = cr_k2_trajectory_;
        Vector<T> pcr_k1                 = pcr_k1_trajectory_;
        Vector<T> pcr_k2                 = pcr_k2_trajectory_;
        Kokkos::View<int*> level_offsets = cr_level_offsets_;
        Vector<T> sub_diagonal_in        = this->sub_diagonal_; // cyclic corner element, unpadded input
        Vector<T> gamma                  = gamma_;

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;

        if (!is_cyclic) {
            const std::size_t scratch_bytes =
                (2ull * static_cast<std::size_t>(n_padded) + 2ull * static_cast<std::size_t>(m)) * sizeof(T);

            TeamPolicy policy(effective_batch_count, n_padded);
            policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

            Kokkos::parallel_for(
                "SolveNonCyclicCRPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                    const int k             = team_member.league_rank();
                    const int batch_idx     = batch_stride * k + batch_offset;
                    const int rank          = team_member.team_rank();
                    const int team_size     = team_member.team_size();
                    const int in_offset     = batch_idx * matrix_dimension;
                    const int padded_offset = batch_idx * n_padded;

                    T* scratch  = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                    T* d[2]     = {scratch, scratch + n_padded};
                    T* pd       = scratch + 2 * n_padded; // length 2*m, ping-pong for PCR sub-solve
                    T* pd0      = pd;
                    T* pd1      = pd + m;
                    T* pdbuf[2] = {pd0, pd1};

                    int cur = 0;
                    for (int i = rank; i < n_padded; i += team_size) {
                        d[cur][i] = (i < matrix_dimension) ? rhs(in_offset + i) : T(0);
                    }
                    team_member.team_barrier();

                    // --- CR forward reduction on d, replaying the stored trajectory. ---
                    for (int L = 0; L < num_cr_levels; ++L) {
                        const int stride       = 1 << (L + 1);
                        const int active_count = n_padded / stride;
                        const int nxt          = 1 - cur;

                        for (int t = rank; t < active_count; t += team_size) {
                            const int i     = stride * t + stride - 1;
                            const int iLeft = i - stride / 2;
                            int iRight      = i + stride / 2;
                            if (iRight > n_padded - 1) {
                                iRight = n_padded - 1;
                            }

                            const T k1_val =
                                cr_k1(static_cast<std::size_t>(batch_idx) * cr_traj_len + level_offsets(L) + t);
                            const T k2_val =
                                cr_k2(static_cast<std::size_t>(batch_idx) * cr_traj_len + level_offsets(L) + t);

                            d[nxt][i] = d[cur][i] - d[cur][iLeft] * k1_val - d[cur][iRight] * k2_val;
                        }
                        for (int i = rank; i < n_padded; i += team_size) {
                            if (!cr_is_active_at_level(i, stride, active_count)) {
                                d[nxt][i] = d[cur][i];
                            }
                        }
                        team_member.team_barrier();
                        cur = nxt;
                    }

                    // --- PCR sub-solve on the m surviving positions. ---
                    const int stride_final = n_padded / m;
                    for (int t = rank; t < m; t += team_size) {
                        const int i = stride_final * t + stride_final - 1;
                        pdbuf[0][t] = d[cur][i];
                    }
                    team_member.team_barrier();

                    int pcur = 0;
                    for (int step = 0; step < num_pcr_steps; ++step) {
                        const int delta = 1 << step;
                        const int pnxt  = 1 - pcur;
                        for (int t = rank; t < m; t += team_size) {
                            int tLeft, tRight;
                            pcr_neighbors(t, delta, m, tLeft, tRight);

                            const T k1_val = pcr_k1(static_cast<std::size_t>(batch_idx) * num_pcr_steps * m +
                                                    static_cast<std::size_t>(step) * m + t);
                            const T k2_val = pcr_k2(static_cast<std::size_t>(batch_idx) * num_pcr_steps * m +
                                                    static_cast<std::size_t>(step) * m + t);

                            pdbuf[pnxt][t] =
                                pdbuf[pcur][t] - pdbuf[pcur][tLeft] * k1_val - pdbuf[pcur][tRight] * k2_val;
                        }
                        team_member.team_barrier();
                        pcur = pnxt;
                    }

                    // Divide by the fully-PCR-reduced diagonal to get x directly at the m
                    // positions, scattering back into the working array (which now holds a mix
                    // of "solved x" at these m positions and "last-forward-touch d" elsewhere,
                    // exactly what CR backward substitution below needs).
                    for (int t = rank; t < m; t += team_size) {
                        const int i = stride_final * t + stride_final - 1;
                        d[cur][i]   = pdbuf[pcur][t] / frozen_b(padded_offset + i);
                    }
                    team_member.team_barrier();

                    // --- CR backward substitution. ---
                    // See the class-level comment for the derivation: neighbors are at
                    // i -+ stride (the FULL stride of the level where i was last updated), not
                    // i -+ stride/2, and each position is solved exactly once, at its true last
                    // forward-reduction touch. The loop runs one level further than the forward
                    // levels (down to a conceptual stride-1 "level -1") to resolve positions
                    // that were never touched by CR forward reduction at all. Positions active
                    // at the coarsest CR level (the m PCR-survivor positions) need no CR
                    // backward step -- PCR's full reduction already solved them above.
                    //
                    // At level L, a position i satisfying level L's active-position formula was
                    // ALSO active at level L+1 -- i.e. already resolved, either by a coarser
                    // level's backward step or directly by PCR -- exactly when its slot index t
                    // is odd (since i = stride*t + stride - 1 also matches level L+1's formula,
                    // with stride doubled, iff t is odd). So only even t need solving here.
                    for (int L = num_cr_levels - 2; L >= -1; --L) {
                        const int stride       = 1 << (L + 1);
                        const int active_count = n_padded / stride;

                        for (int t = rank; t < active_count; t += team_size) {
                            if ((t & 1) != 0) {
                                continue;
                            }
                            const int i = stride * t + stride - 1;
                            int iLeft   = i - stride;
                            if (iLeft < 0) {
                                iLeft = 0;
                            }
                            int iRight = i + stride;
                            if (iRight > n_padded - 1) {
                                iRight = n_padded - 1;
                            }

                            const T rhs_val = d[cur][i];
                            const T xLeft   = d[cur][iLeft];
                            const T xRight  = d[cur][iRight];
                            d[cur][i] =
                                (rhs_val - frozen_a(padded_offset + i) * xLeft - frozen_c(padded_offset + i) * xRight) /
                                frozen_b(padded_offset + i);
                        }
                        team_member.team_barrier();
                    }

                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        rhs(in_offset + i) = d[cur][i];
                    }
                });
        }
        else {
            // Cyclic solve: run the same forward/PCR/backward pipeline twice per kernel launch
            // against the same stored trajectories -- once for the rhs-derived vector, once for
            // the Sherman-Morrison auxiliary vector (gamma at position 0, cyclic_corner at
            // position matrix_dimension - 1, zero elsewhere including the padding region) --
            // then combine via the same rank-1 (Sherman-Morrison) formula used by the existing
            // solvers. This logic is unchanged from BatchedTridiagonalSolverPCR/Thomas; only the
            // underlying "solve a plain tridiagonal system" primitive has changed.
            const std::size_t scratch_bytes =
                (4ull * static_cast<std::size_t>(n_padded) + 4ull * static_cast<std::size_t>(m)) * sizeof(T);

            TeamPolicy policy(effective_batch_count, n_padded);
            policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

            Kokkos::parallel_for(
                "SolveCyclicCRPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                    const int k             = team_member.league_rank();
                    const int batch_idx     = batch_stride * k + batch_offset;
                    const int rank          = team_member.team_rank();
                    const int team_size     = team_member.team_size();
                    const int in_offset     = batch_idx * matrix_dimension;
                    const int padded_offset = batch_idx * n_padded;

                    T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                    // Layout: [d_rhs(0) | d_buf(0) | d_rhs(1) | d_buf(1) | pd_rhs(0) | pd_buf(0)
                    //          | pd_rhs(1) | pd_buf(1)].
                    T* d_rhs[2]  = {scratch, scratch + 2 * n_padded};
                    T* d_buf[2]  = {scratch + n_padded, scratch + 3 * n_padded};
                    T* p_base    = scratch + 4 * n_padded;
                    T* pd_rhs[2] = {p_base, p_base + 2 * m};
                    T* pd_buf[2] = {p_base + m, p_base + 3 * m};

                    const T cyclic_corner_element = sub_diagonal_in(in_offset + matrix_dimension - 1);

                    int cur = 0;
                    for (int i = rank; i < n_padded; i += team_size) {
                        if (i < matrix_dimension) {
                            d_rhs[cur][i] = rhs(in_offset + i);
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
                    team_member.team_barrier();

                    for (int L = 0; L < num_cr_levels; ++L) {
                        const int stride       = 1 << (L + 1);
                        const int active_count = n_padded / stride;
                        const int nxt          = 1 - cur;

                        for (int t = rank; t < active_count; t += team_size) {
                            const int i     = stride * t + stride - 1;
                            const int iLeft = i - stride / 2;
                            int iRight      = i + stride / 2;
                            if (iRight > n_padded - 1) {
                                iRight = n_padded - 1;
                            }

                            const T k1_val =
                                cr_k1(static_cast<std::size_t>(batch_idx) * cr_traj_len + level_offsets(L) + t);
                            const T k2_val =
                                cr_k2(static_cast<std::size_t>(batch_idx) * cr_traj_len + level_offsets(L) + t);

                            d_rhs[nxt][i] = d_rhs[cur][i] - d_rhs[cur][iLeft] * k1_val - d_rhs[cur][iRight] * k2_val;
                            d_buf[nxt][i] = d_buf[cur][i] - d_buf[cur][iLeft] * k1_val - d_buf[cur][iRight] * k2_val;
                        }
                        for (int i = rank; i < n_padded; i += team_size) {
                            if (!cr_is_active_at_level(i, stride, active_count)) {
                                d_rhs[nxt][i] = d_rhs[cur][i];
                                d_buf[nxt][i] = d_buf[cur][i];
                            }
                        }
                        team_member.team_barrier();
                        cur = nxt;
                    }

                    const int stride_final = n_padded / m;
                    for (int t = rank; t < m; t += team_size) {
                        const int i  = stride_final * t + stride_final - 1;
                        pd_rhs[0][t] = d_rhs[cur][i];
                        pd_buf[0][t] = d_buf[cur][i];
                    }
                    team_member.team_barrier();

                    int pcur = 0;
                    for (int step = 0; step < num_pcr_steps; ++step) {
                        const int delta = 1 << step;
                        const int pnxt  = 1 - pcur;
                        for (int t = rank; t < m; t += team_size) {
                            int tLeft, tRight;
                            pcr_neighbors(t, delta, m, tLeft, tRight);

                            const T k1_val = pcr_k1(static_cast<std::size_t>(batch_idx) * num_pcr_steps * m +
                                                    static_cast<std::size_t>(step) * m + t);
                            const T k2_val = pcr_k2(static_cast<std::size_t>(batch_idx) * num_pcr_steps * m +
                                                    static_cast<std::size_t>(step) * m + t);

                            pd_rhs[pnxt][t] =
                                pd_rhs[pcur][t] - pd_rhs[pcur][tLeft] * k1_val - pd_rhs[pcur][tRight] * k2_val;
                            pd_buf[pnxt][t] =
                                pd_buf[pcur][t] - pd_buf[pcur][tLeft] * k1_val - pd_buf[pcur][tRight] * k2_val;
                        }
                        team_member.team_barrier();
                        pcur = pnxt;
                    }

                    for (int t = rank; t < m; t += team_size) {
                        const int i   = stride_final * t + stride_final - 1;
                        d_rhs[cur][i] = pd_rhs[pcur][t] / frozen_b(padded_offset + i);
                        d_buf[cur][i] = pd_buf[pcur][t] / frozen_b(padded_offset + i);
                    }
                    team_member.team_barrier();

                    for (int L = num_cr_levels - 2; L >= -1; --L) {
                        const int stride       = 1 << (L + 1);
                        const int active_count = n_padded / stride;

                        for (int t = rank; t < active_count; t += team_size) {
                            if ((t & 1) != 0) {
                                continue;
                            }
                            const int i = stride * t + stride - 1;
                            int iLeft   = i - stride;
                            if (iLeft < 0) {
                                iLeft = 0;
                            }
                            int iRight = i + stride;
                            if (iRight > n_padded - 1) {
                                iRight = n_padded - 1;
                            }

                            const T fa = frozen_a(padded_offset + i);
                            const T fb = frozen_b(padded_offset + i);
                            const T fc = frozen_c(padded_offset + i);

                            const T rhs_val_rhs = d_rhs[cur][i];
                            const T rhs_val_buf = d_buf[cur][i];

                            d_rhs[cur][i] = (rhs_val_rhs - fa * d_rhs[cur][iLeft] - fc * d_rhs[cur][iRight]) / fb;
                            d_buf[cur][i] = (rhs_val_buf - fa * d_buf[cur][iLeft] - fc * d_buf[cur][iRight]) / fb;
                        }
                        team_member.team_barrier();
                    }

                    // Sherman-Morrison reconstruction, identical to the existing PCR/Thomas
                    // implementations.
                    const T dot_product_x_v =
                        d_rhs[cur][0] + cyclic_corner_element / gamma(batch_idx) * d_rhs[cur][matrix_dimension - 1];
                    const T dot_product_u_v =
                        d_buf[cur][0] + cyclic_corner_element / gamma(batch_idx) * d_buf[cur][matrix_dimension - 1];
                    const T factor = dot_product_x_v / (T(1) + dot_product_u_v);

                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        rhs(in_offset + i) = d_rhs[cur][i] - factor * d_buf[cur][i];
                    }
                });
        }
        Kokkos::fence();
    }

    /**
     * @brief Solves systems whose matrix has already been reduced to diagonal form.
     *
     * Unrelated to the CR/PCR reduction machinery -- unchanged in behavior from the existing
     * implementations, operating on the original (unpadded) main_diagonal_ input plus the
     * cyclic gamma correction.
     */
    void solve_diagonal(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1) override
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        const int effective_batch_count = (this->batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        const int matrix_dimension = this->matrix_dimension_;
        const bool is_cyclic       = this->is_cyclic_;
        Vector<T> main_diagonal    = this->main_diagonal_;
        Vector<T> gamma            = gamma_;

        using MDPolicy = Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>;
        MDPolicy policy({0, 0}, {effective_batch_count, matrix_dimension});

        Kokkos::parallel_for(
            "SolveDiagonalCRPCR", policy, KOKKOS_LAMBDA(const int k, const int i) {
                const int batch_idx = batch_stride * k + batch_offset;
                const int offset    = batch_idx * matrix_dimension;

                if (is_cyclic && i == 0) {
                    rhs(offset) /= main_diagonal(offset) + gamma(batch_idx);
                }
                else {
                    rhs(offset + i) /= main_diagonal(offset + i);
                }
            });

        Kokkos::fence();
    }

    // Exposed for tests (e.g. the memory-footprint regression check).
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

private:
    int matrix_dimension_padded_; // next_power_of_two(matrix_dimension_)
    int intermediate_system_size_; // m
    int num_cr_levels_; // log2(matrix_dimension_padded_ / m)
    int num_pcr_steps_; // log2(m)
    int cr_trajectory_length_per_batch_; // ~ matrix_dimension_padded_

    // Frozen forward-reduction state (O(n) each); see class-level comment.
    AllocatableVector<T> frozen_sub_lower_; // "a"
    AllocatableVector<T> frozen_main_diagonal_; // "b"
    AllocatableVector<T> frozen_sub_upper_; // "c"

    // Compact CR multiplier trajectory (O(n) each) and per-level prefix-sum offsets.
    AllocatableVector<T> cr_k1_trajectory_;
    AllocatableVector<T> cr_k2_trajectory_;
    Kokkos::View<int*> cr_level_offsets_;

    // PCR trajectory for the small intermediate system (O(m log m) each).
    AllocatableVector<T> pcr_k1_trajectory_;
    AllocatableVector<T> pcr_k2_trajectory_;

    Vector<T> gamma_;
    bool is_factorized_;
};

} // namespace gmgpolar