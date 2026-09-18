#pragma once

/**
 * @brief Batched tridiagonal solver using hybrid Cyclic Reduction (CR) +
 *        Parallel Cyclic Reduction (PCR).
 *
 * This solver is built directly on top of BatchedTridiagonalSolverCR's
 * in-place CR factorization/backward-substitution algebra and storage
 * model, and on top of BatchedTridiagonalSolverPCR's compact reduction
 * recurrence. It is not a redesign of either: it runs the *existing* CR
 * levels while the active system is large, and switches to the *existing*
 * PCR recurrence -- restricted to a small compact survivor system -- once
 * the number of surviving equations is small.
 *
 * -------------------------------------------------------------------------
 * Algorithm
 * -------------------------------------------------------------------------
 *
 *   original system (size n)
 *        |
 *        v
 *   CR forward reduction, in place, exactly as in BatchedTridiagonalSolverCR
 *        |   (stops after L levels, as soon as the number of CR survivors
 *        |    m satisfies m <= PCR_TARGET_SIZE)
 *        v
 *   compact survivor system (size m <= 128)
 *        |
 *        v
 *   compact PCR factorization/solve (team scratch, O(m) per buffer)
 *        |
 *        v
 *   survivor solutions scattered back into their original CR positions
 *        |
 *        v
 *   CR backward substitution, in place, exactly as in BatchedTridiagonalSolverCR
 *        |
 *        v
 *   full solution
 */

#include <Kokkos_Core.hpp>
#include <cstddef>
#include <stdexcept>

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

namespace gmgpolar
{

template <typename T>
class BatchedTridiagonalSolverCRPCR : public BatchedTridiagonalSolverBase<T>
{
public:
    // Fixed PCR target: CR runs until the survivor count is <= this value.
    // No runtime cutoff parameter exists; this constant is isolated here
    // so a future version can change it without touching the CR logic.
    static constexpr int PCR_TARGET_SIZE = 128;

    BatchedTridiagonalSolverCRPCR(int matrix_dimension, int batch_count, bool is_cyclic = true)
        : BatchedTridiagonalSolverBase<T>(matrix_dimension, batch_count, is_cyclic)
        , num_cr_levels_(compute_num_cr_levels(matrix_dimension))
        , survivor_base_(num_cr_levels_ > 0 ? ((1 << num_cr_levels_) - 1) : 0)
        , survivor_stride_(num_cr_levels_ > 0 ? (1 << num_cr_levels_) : 1)
        , m_(compute_survivor_count(matrix_dimension, num_cr_levels_))
        , pcr_num_steps_(compute_pcr_steps(m_))
        , q_right_trajectory_("BatchedTridiagonalSolverCRPCR::q_right_trajectory",
                              static_cast<std::size_t>(matrix_dimension) * static_cast<std::size_t>(batch_count))
        , cyclic_corner_("BatchedTridiagonalSolverCRPCR::cyclic_corner", is_cyclic ? batch_count : 0)
        , gamma_("BatchedTridiagonalSolverCRPCR::gamma", is_cyclic ? batch_count : 0)
        , pcr_k1_trajectory_("BatchedTridiagonalSolverCRPCR::pcr_k1_trajectory",
                             static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(m_) *
                                 static_cast<std::size_t>(pcr_num_steps_))
        , pcr_k2_trajectory_("BatchedTridiagonalSolverCRPCR::pcr_k2_trajectory",
                             static_cast<std::size_t>(batch_count) * static_cast<std::size_t>(m_) *
                                 static_cast<std::size_t>(pcr_num_steps_))
        , is_factorized_(false)
    {
        assign(q_right_trajectory_, T(0));
        assign(cyclic_corner_, T(0));
        assign(gamma_, T(0));
        assign(pcr_k1_trajectory_, T(0));
        assign(pcr_k2_trajectory_, T(0));
    }

    // Number of CR levels actually executed (0 when matrix_dimension <= PCR_TARGET_SIZE).
    int crLevels() const
    {
        return num_cr_levels_;
    }

    // Number of CR survivors / compact PCR dimension (<= PCR_TARGET_SIZE).
    int survivorCount() const
    {
        return m_;
    }

    // Number of compact PCR reduction steps needed to reduce m equations, 0 if m<=1.
    int pcrSteps() const
    {
        return pcr_num_steps_;
    }

    /* -------------------------------------------------------------- */
    /* Setup: CR forward factorization (existing algebra) down to the  */
    /* survivor boundary, then compact PCR factorization of the m<=128 */
    /* survivor system, writing the final diagonal back into the       */
    /* survivor slots of main_diagonal_.                                */
    /* -------------------------------------------------------------- */
    void setup() override
    {
        const int n               = this->matrix_dimension_;
        const int L               = num_cr_levels_;
        const int survivor_base   = survivor_base_;
        const int survivor_stride = survivor_stride_;
        const int m               = m_;
        const int pcr_steps       = pcr_num_steps_;
        const bool is_cyclic      = this->is_cyclic_;

        Vector<T> main_diagonal      = this->main_diagonal_;
        Vector<T> sub_diagonal       = this->sub_diagonal_;
        Vector<T> q_right_trajectory = q_right_trajectory_;
        Vector<T> gamma              = gamma_;
        Vector<T> cyclic_corner      = cyclic_corner_;
        Vector<T> pcr_k1             = pcr_k1_trajectory_;
        Vector<T> pcr_k2             = pcr_k2_trajectory_;

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;

        // n == 1: the cyclic-corner storage slot IS index 0
        // (matrix_dimension - 1 == 0), so the generic Sherman-Morrison
        // derivation is not meaningful. BatchedTridiagonalSolverPCR
        // explicitly special-cases matrix_dimension == 1 the same way:
        // it is treated as a plain 1x1 diagonal system regardless of
        // is_cyclic. CRPCR matches that established convention exactly
        // (verified against a dense reference in the standalone harness).
        if (n == 1) {
            is_factorized_ = true;
            return;
        }

        const std::size_t scratch_bytes = 4ull * static_cast<std::size_t>(m) * sizeof(T);

        TeamPolicy policy(this->batch_count_, Kokkos::AUTO);
        policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));

        Kokkos::parallel_for(
            "SetupCRPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                const int batch_idx = team_member.league_rank();
                const int offset    = batch_idx * n;
                const int team_size = team_member.team_size();
                const int rank      = team_member.team_rank();

                // -------------------------------------------------- //
                // Sherman-Morrison Adjustment (cyclic systems only), //
                // identical to BatchedTridiagonalSolverCR::setup().  //
                // Saved BEFORE CR can overwrite sub_diagonal_.       //
                // -------------------------------------------------- //
                if (is_cyclic) {
                    if (rank == 0) {
                        const T corner_element   = sub_diagonal(offset + n - 1);
                        cyclic_corner(batch_idx) = corner_element;
                        gamma(batch_idx)         = -main_diagonal(offset + 0);
                        main_diagonal(offset + 0) -= gamma(batch_idx);
                        main_diagonal(offset + n - 1) -= corner_element * corner_element / gamma(batch_idx);
                    }
                    team_member.team_barrier();
                }

                // --------------------------------------------------- //
                // CR forward levels 0 .. L-1: IDENTICAL arithmetic to //
                // BatchedTridiagonalSolverCR::setup(), just bounded   //
                // to L levels instead of running to a single root.    //
                // --------------------------------------------------- //
                for (int step = 0; step < L; ++step) {
                    const int d = 1 << step;

                    // Phase A: eliminate i = (d-1) + 2*d*k.
                    for (int k = rank;; k += team_size) {
                        const int i = (d - 1) + 2 * d * k;
                        if (i >= n) {
                            break;
                        }

                        const bool has_left  = (i - d) >= 0;
                        const bool has_right = (i + d) < n;

                        const T b_i = main_diagonal(offset + i);
                        const T a_i = has_left ? sub_diagonal(offset + i - d) : T(0);
                        const T c_i = has_right ? sub_diagonal(offset + i) : T(0);

                        const T inv_b  = T(1) / b_i;
                        const T qLeft  = a_i / b_i;
                        const T qRight = c_i / b_i;

                        main_diagonal(offset + i)      = inv_b;
                        sub_diagonal(offset + i)       = qLeft;
                        q_right_trajectory(offset + i) = qRight;
                    }

                    team_member.team_barrier();

                    // Phase B: survivors s = (2*d-1) + 2*d*k.
                    for (int k = rank;; k += team_size) {
                        const int s = (2 * d - 1) + 2 * d * k;
                        if (s >= n) {
                            break;
                        }

                        const int eL      = s - d;
                        const int eR      = s + d;
                        const bool has_eR = eR < n;

                        const T qRightEL = q_right_trajectory(offset + eL);

                        const T c_s_raw = sub_diagonal(offset + s);
                        const T b_s_raw = main_diagonal(offset + s);

                        const T term1 = qRightEL * qRightEL / main_diagonal(offset + eL);

                        T term2 = T(0);
                        T c_new = T(0);
                        if (has_eR) {
                            const T qLeftER  = sub_diagonal(offset + eR);
                            const T qRightER = q_right_trajectory(offset + eR);
                            term2            = qLeftER * c_s_raw;
                            c_new            = -c_s_raw * qRightER;
                        }

                        main_diagonal(offset + s) = b_s_raw - term1 - term2;
                        sub_diagonal(offset + s)  = c_new;
                    }

                    team_member.team_barrier();
                }

                // ------------------------------------------------- //
                // Compact gather: survivors r_j = survivor_base +   //
                // j*survivor_stride hold "current reduced diagonal" //
                // (main_diagonal) and "current reduced right        //
                // coefficient" (sub_diagonal), per the persistent   //
                // storage contract.                                 //
                // ------------------------------------------------- //
                T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                // Layout: [e0 | b0 | e1 | b1], each length m.
                T* e[2] = {scratch, scratch + 2 * m};
                T* b[2] = {scratch + m, scratch + 3 * m};

                for (int j = rank; j < m; j += team_size) {
                    const int orig = survivor_base + j * survivor_stride;
                    b[0][j]        = main_diagonal(offset + orig);
                    e[0][j]        = sub_diagonal(offset + orig);
                }
                team_member.team_barrier();

                // When L==0 (n<=PCR_TARGET_SIZE: no CR ran), the survivor
                // set is the raw original system and the last equation's
                // sub_diagonal_ slot is either unused (non-cyclic) or holds
                // the cyclic corner (cyclic) -- never a real right-neighbor
                // coupling. Zero it explicitly, exactly as
                // BatchedTridiagonalSolverPCR::setup() already does
                // (e[cur][i] = (i==n-1) ? 0 : sub_diagonal(...)). When
                // L>=1 this is already guaranteed by CR's own Phase B
                // (c_new==0 whenever the survivor has no right CR
                // neighbor, which is exactly true for the last survivor),
                // so no separate action is required in that case.
                if (L == 0 && rank == 0 && m > 0) {
                    e[0][m - 1] = T(0);
                }
                team_member.team_barrier();

                // ----------------------------------------------------- //
                // Compact PCR factorization: IDENTICAL recurrence to    //
                // BatchedTridiagonalSolverPCR::setup(), re-scoped from  //
                // size n to size m. Symmetry is used: only e[] (the     //
                // right coefficient) is stored per level; the left      //
                // coefficient a(j) = c(j-delta) = e[cur][j-delta] is    //
                // reconstructed, never stored separately.               //
                // Two multiplier trajectories (k1, k2) are stored       //
                // because k1(j) != k2(j-delta) in general (Section 6):  //
                // this mirrors the already-proven-correct convention of //
                // the reference PCR solver rather than inventing a new, //
                // unproven single-array compression.                    //
                // ----------------------------------------------------- //
                int cur = 0;
                for (int step = 0; step < pcr_steps; ++step) {
                    const int delta = 1 << step;

                    for (int j = rank; j < m; j += team_size) {
                        const bool has_left  = j >= delta;
                        const bool has_right = j + delta < m;
                        const int left       = j - delta;
                        const int right      = j + delta;

                        // The compact system is symmetric. Therefore the left
                        // coefficient a(j) is c(j-1) and, at this PCR level,
                        // a(j) = c(j-delta). Boundary terms must be zero rather
                        // than clamped to an existing equation.
                        const T a_j     = has_left ? e[cur][left] : T(0);
                        const T c_j     = has_right ? e[cur][j] : T(0);
                        const T c_right = has_right ? e[cur][right] : T(0);

                        const T k1_val = has_left ? a_j / b[cur][left] : T(0);
                        const T k2_val = has_right ? c_j / b[cur][right] : T(0);

                        const std::size_t traj = static_cast<std::size_t>(batch_idx) *
                                                     static_cast<std::size_t>(pcr_steps) * static_cast<std::size_t>(m) +
                                                 static_cast<std::size_t>(step) * static_cast<std::size_t>(m) +
                                                 static_cast<std::size_t>(j);

                        pcr_k1(traj) = k1_val;
                        pcr_k2(traj) = k2_val;

                        const int nxt = 1 - cur;
                        e[nxt][j]     = has_right ? -c_right * k2_val : T(0);
                        b[nxt][j]     = b[cur][j] - a_j * k1_val - c_j * k2_val;
                    }

                    team_member.team_barrier();
                    cur = 1 - cur;
                }

                // ---------------------------------------------------- //
                // Write the final PCR diagonal back into the survivor  //
                // positions of main_diagonal_ (Section 10). No         //
                // persistent compact pcr_main_diagonal_ array is       //
                // created. The final PCR off-diagonal is not needed    //
                // by solve() and is not persisted (Section 11).        //
                // Eliminated (non-survivor) main_diagonal_ entries are //
                // untouched by this loop.                              //
                // ---------------------------------------------------- //
                for (int j = rank; j < m; j += team_size) {
                    const int orig               = survivor_base + j * survivor_stride;
                    main_diagonal(offset + orig) = b[cur][j];
                }
            });

        Kokkos::fence();
        is_factorized_ = true;
    }

    /* --------------------------------------------------------------- */
    /* Solve: CR forward RHS reduction (existing algebra) -> compact   */
    /* survivor gather -> compact PCR solve (replay of the stored      */
    /* factorization) -> scatter -> CR backward substitution (existing */
    /* algebra, unchanged).                                            */
    /* --------------------------------------------------------------- */
    void solve(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1) override
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }
        if (batch_stride <= 0) {
            throw std::invalid_argument("Error: batch_stride must be positive.");
        }
        if (batch_offset < 0 || batch_offset > this->batch_count_) {
            throw std::invalid_argument("Error: batch_offset out of range.");
        }

        const int effective_batch_count = compute_effective_batch_count(batch_offset, batch_stride);
        if (effective_batch_count == 0) {
            return;
        }

        const int n               = this->matrix_dimension_;
        const int L               = num_cr_levels_;
        const int survivor_base   = survivor_base_;
        const int survivor_stride = survivor_stride_;
        const int m               = m_;
        const int pcr_steps       = pcr_num_steps_;
        const bool is_cyclic      = this->is_cyclic_;

        Vector<T> main_diagonal      = this->main_diagonal_;
        Vector<T> sub_diagonal       = this->sub_diagonal_;
        Vector<T> q_right_trajectory = q_right_trajectory_;
        Vector<T> gamma              = gamma_;
        Vector<T> cyclic_corner      = cyclic_corner_;
        Vector<T> pcr_k1             = pcr_k1_trajectory_;
        Vector<T> pcr_k2             = pcr_k2_trajectory_;

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;

        if (n == 1) {
            // Matches BatchedTridiagonalSolverPCR's n==1 special case:
            // plain diagonal solve, regardless of is_cyclic.
            Kokkos::parallel_for(
                "SolveCRPCRTrivial", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, effective_batch_count),
                KOKKOS_LAMBDA(const int k) {
                    const int batch_idx = batch_stride * k + batch_offset;
                    rhs(batch_idx) /= main_diagonal(batch_idx);
                });
            Kokkos::fence();
            return;
        }

        if (!is_cyclic) {
            const std::size_t scratch_bytes = 2ull * static_cast<std::size_t>(m) * sizeof(T);

            TeamPolicy policy(effective_batch_count, Kokkos::AUTO);
            policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));

            Kokkos::parallel_for(
                "SolveCRPCRNonCyclic", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                    const int k         = team_member.league_rank();
                    const int batch_idx = batch_stride * k + batch_offset;
                    const int offset    = batch_idx * n;
                    const int team_size = team_member.team_size();
                    const int rank      = team_member.team_rank();

                    // -------------------- //
                    // CR Forward RHS Reduction (levels 0..L-1), identical //
                    // algebra to BatchedTridiagonalSolverCR::solve().     //
                    // -------------------- //
                    for (int step = 0; step < L; ++step) {
                        const int d = 1 << step;

                        for (int kk = rank;; kk += team_size) {
                            const int s = (2 * d - 1) + 2 * d * kk;
                            if (s >= n) {
                                break;
                            }

                            const int eL      = s - d;
                            const int eR      = s + d;
                            const bool has_eR = eR < n;

                            const T factor_left  = q_right_trajectory(offset + eL);
                            const T factor_right = has_eR ? sub_diagonal(offset + eR) : T(0);
                            const T rhs_eR       = has_eR ? rhs(offset + eR) : T(0);

                            rhs(offset + s) -= factor_left * rhs(offset + eL) + factor_right * rhs_eR;
                        }

                        team_member.team_barrier();
                    }

                    // ------------- //
                    // Compact survivor gather //
                    // ------------- //
                    T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                    T* d[2]    = {scratch, scratch + m};

                    int cur = 0;
                    for (int j = rank; j < m; j += team_size) {
                        const int orig = survivor_base + j * survivor_stride;
                        d[cur][j]      = rhs(offset + orig);
                    }
                    team_member.team_barrier();

                    // ----------------------------------------------------- //
                    // Compact PCR solve (trajectory replay), identical to   //
                    // BatchedTridiagonalSolverPCR::solve(), re-scoped to m. //
                    // ----------------------------------------------------- //
                    for (int step = 0; step < pcr_steps; ++step) {
                        const int delta = 1 << step;

                        for (int j = rank; j < m; j += team_size) {
                            const bool has_left  = j >= delta;
                            const bool has_right = j + delta < m;
                            const int left       = j - delta;
                            const int right      = j + delta;

                            const std::size_t traj = static_cast<std::size_t>(batch_idx) *
                                                         static_cast<std::size_t>(pcr_steps) *
                                                         static_cast<std::size_t>(m) +
                                                     static_cast<std::size_t>(step) * static_cast<std::size_t>(m) +
                                                     static_cast<std::size_t>(j);

                            const T k1_val = pcr_k1(traj);
                            const T k2_val = pcr_k2(traj);

                            const int nxt = 1 - cur;

                            d[nxt][j] = d[cur][j] - (has_left ? k1_val * d[cur][left] : T(0)) -
                                        (has_right ? k2_val * d[cur][right] : T(0));
                        }

                        team_member.team_barrier();
                        cur = 1 - cur;
                    }

                    // PCR is a complete solve: divide by the final diagonal
                    // and scatter directly into the original survivor
                    // positions. No PCR backward-substitution phase.
                    for (int j = rank; j < m; j += team_size) {
                        const int orig     = survivor_base + j * survivor_stride;
                        rhs(offset + orig) = d[cur][j] / main_diagonal(offset + orig);
                    }
                    team_member.team_barrier();

                    // --------------------------------------------------- //
                    // CR Backward Substitution (levels L-1..0), identical //
                    // algebra to BatchedTridiagonalSolverCR::solve().     //
                    // --------------------------------------------------- //
                    for (int step = L - 1; step >= 0; --step) {
                        const int d = 1 << step;

                        for (int kk = rank;; kk += team_size) {
                            const int i = (d - 1) + 2 * d * kk;
                            if (i >= n) {
                                break;
                            }

                            const bool has_left  = (i - d) >= 0;
                            const bool has_right = (i + d) < n;

                            const T x_left  = has_left ? rhs(offset + i - d) : T(0);
                            const T x_right = has_right ? rhs(offset + i + d) : T(0);

                            rhs(offset + i) = rhs(offset + i) * main_diagonal(offset + i) -
                                              sub_diagonal(offset + i) * x_left -
                                              q_right_trajectory(offset + i) * x_right;
                        }

                        team_member.team_barrier();
                    }
                });
        }
        else {
            // Cyclic solve: the actual RHS and the Sherman-Morrison
            // auxiliary vector are carried through the identical
            // CR-forward -> compact-PCR -> CR-backward pipeline together,
            // reusing the single stored factorization, then combined via
            // the existing Sherman-Morrison reconstruction. The auxiliary
            // vector uses O(n) team scratch (the same pattern already used
            // by BatchedTridiagonalSolverCR::solve()'s cyclic branch)
            // rather than a permanent batch_count*n array.
            const std::size_t aux_bytes     = static_cast<std::size_t>(n) * sizeof(T);
            const std::size_t compact_bytes = 4ull * static_cast<std::size_t>(m) * sizeof(T);
            const std::size_t scratch_bytes = aux_bytes + compact_bytes;

            TeamPolicy policy(effective_batch_count, Kokkos::AUTO);
            policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));

            Kokkos::parallel_for(
                "SolveCRPCRCyclic", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                    const int k         = team_member.league_rank();
                    const int batch_idx = batch_stride * k + batch_offset;
                    const int offset    = batch_idx * n;
                    const int team_size = team_member.team_size();
                    const int rank      = team_member.team_rank();

                    T* raw      = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                    T* buffer   = raw; // O(n) Sherman-Morrison auxiliary vector.
                    T* compact  = raw + n; // O(4m) compact PCR scratch.
                    T* d_rhs[2] = {compact, compact + 2 * m};
                    T* d_buf[2] = {compact + m, compact + 3 * m};

                    // ------------------------------------------------- //
                    // Initialize the Sherman-Morrison auxiliary vector: //
                    // u = gamma * e_0 + corner * e_{n-1}, identical to  //
                    // BatchedTridiagonalSolverCR::solve().              //
                    // ------------------------------------------------- //
                    for (int i = rank; i < n; i += team_size) {
                        buffer[i] = T(0);
                    }
                    team_member.team_barrier();
                    if (rank == 0) {
                        buffer[0] = gamma(batch_idx);
                        if (n > 1) {
                            buffer[n - 1] = cyclic_corner(batch_idx);
                        }
                    }
                    team_member.team_barrier();

                    // ------------------------------------------------- //
                    // CR Forward RHS Reduction, actual RHS + auxiliary, //
                    // identical algebra to BatchedTridiagonalSolverCR.  //
                    // ------------------------------------------------- //
                    for (int step = 0; step < L; ++step) {
                        const int d = 1 << step;

                        for (int kk = rank;; kk += team_size) {
                            const int s = (2 * d - 1) + 2 * d * kk;
                            if (s >= n) {
                                break;
                            }

                            const int eL      = s - d;
                            const int eR      = s + d;
                            const bool has_eR = eR < n;

                            const T factor_left  = q_right_trajectory(offset + eL);
                            const T factor_right = has_eR ? sub_diagonal(offset + eR) : T(0);
                            const T rhs_eR       = has_eR ? rhs(offset + eR) : T(0);
                            const T buf_eR       = has_eR ? buffer[eR] : T(0);

                            rhs(offset + s) -= factor_left * rhs(offset + eL) + factor_right * rhs_eR;
                            buffer[s] -= factor_left * buffer[eL] + factor_right * buf_eR;
                        }

                        team_member.team_barrier();
                    }

                    // --------------------------------- //
                    // Compact gather (both quantities). //
                    // --------------------------------- //
                    int cur = 0;
                    for (int j = rank; j < m; j += team_size) {
                        const int orig = survivor_base + j * survivor_stride;
                        d_rhs[cur][j]  = rhs(offset + orig);
                        d_buf[cur][j]  = buffer[orig];
                    }
                    team_member.team_barrier();

                    // ----------------------------------------------- //
                    // Compact PCR solve (trajectory replay) for both. //
                    // ----------------------------------------------- //
                    for (int step = 0; step < pcr_steps; ++step) {
                        const int delta = 1 << step;

                        for (int j = rank; j < m; j += team_size) {
                            const bool has_left  = j >= delta;
                            const bool has_right = j + delta < m;
                            const int left       = j - delta;
                            const int right      = j + delta;

                            const std::size_t traj = static_cast<std::size_t>(batch_idx) *
                                                         static_cast<std::size_t>(pcr_steps) *
                                                         static_cast<std::size_t>(m) +
                                                     static_cast<std::size_t>(step) * static_cast<std::size_t>(m) +
                                                     static_cast<std::size_t>(j);

                            const T k1_val = pcr_k1(traj);
                            const T k2_val = pcr_k2(traj);

                            const int nxt = 1 - cur;
                            d_rhs[nxt][j] = d_rhs[cur][j] - (has_left ? k1_val * d_rhs[cur][left] : T(0)) -
                                            (has_right ? k2_val * d_rhs[cur][right] : T(0));
                            d_buf[nxt][j] = d_buf[cur][j] - (has_left ? k1_val * d_buf[cur][left] : T(0)) -
                                            (has_right ? k2_val * d_buf[cur][right] : T(0));
                        }

                        team_member.team_barrier();
                        cur = 1 - cur;
                    }

                    // Scatter both survivor solutions.
                    for (int j = rank; j < m; j += team_size) {
                        const int orig     = survivor_base + j * survivor_stride;
                        rhs(offset + orig) = d_rhs[cur][j] / main_diagonal(offset + orig);
                        buffer[orig]       = d_buf[cur][j] / main_diagonal(offset + orig);
                    }
                    team_member.team_barrier();

                    // ---------------------------------------------------- //
                    // CR Backward Substitution for both, identical algebra //
                    // to BatchedTridiagonalSolverCR::solve().              //
                    // ---------------------------------------------------- //
                    for (int step = L - 1; step >= 0; --step) {
                        const int d = 1 << step;

                        for (int kk = rank;; kk += team_size) {
                            const int i = (d - 1) + 2 * d * kk;
                            if (i >= n) {
                                break;
                            }

                            const bool has_left  = (i - d) >= 0;
                            const bool has_right = (i + d) < n;

                            const T x_left  = has_left ? rhs(offset + i - d) : T(0);
                            const T x_right = has_right ? rhs(offset + i + d) : T(0);
                            const T v_left  = has_left ? buffer[i - d] : T(0);
                            const T v_right = has_right ? buffer[i + d] : T(0);

                            const T invB   = main_diagonal(offset + i);
                            const T qLeft  = sub_diagonal(offset + i);
                            const T qRight = q_right_trajectory(offset + i);

                            rhs(offset + i) = rhs(offset + i) * invB - qLeft * x_left - qRight * x_right;
                            buffer[i]       = buffer[i] * invB - qLeft * v_left - qRight * v_right;
                        }

                        team_member.team_barrier();
                    }

                    // -------------------------------------------------- //
                    // Sherman-Morrison Reconstruction, identical formula //
                    // to BatchedTridiagonalSolverCR::solve().            //
                    // -------------------------------------------------- //
                    const T corner = cyclic_corner(batch_idx);
                    const T g      = gamma(batch_idx);

                    const T dot_product_x_v = rhs(offset + 0) + corner / g * rhs(offset + n - 1);
                    const T dot_product_u_v = buffer[0] + corner / g * buffer[n - 1];
                    const T factor          = dot_product_x_v / (T(1) + dot_product_u_v);

                    for (int i = rank; i < n; i += team_size) {
                        rhs(offset + i) -= factor * buffer[i];
                    }
                });
        }
        Kokkos::fence();
    }

    /* ---------------------------- */
    /* Solve: Diagonal Scaling Only */
    /* ---------------------------- */
    // Does NOT run CR or PCR. Preserves BatchedTridiagonalSolverCR's
    // semantics, generalized from "the single CR root" to "the m CR
    // survivor positions": eliminated (non-survivor) positions store an
    // inverse pivot (1/b_i); survivor positions store the raw final
    // diagonal directly (no PCR backward phase exists to undo, matching
    // BatchedTridiagonalSolverPCR's own convention for its stored root).
    void solve_diagonal(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1) override
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }
        if (batch_stride <= 0) {
            throw std::invalid_argument("Error: batch_stride must be positive.");
        }
        if (batch_offset < 0 || batch_offset > this->batch_count_) {
            throw std::invalid_argument("Error: batch_offset out of range.");
        }

        const int effective_batch_count = compute_effective_batch_count(batch_offset, batch_stride);
        if (effective_batch_count == 0) {
            return;
        }

        const int n               = this->matrix_dimension_;
        const int L               = num_cr_levels_;
        const int survivor_base   = survivor_base_;
        const int survivor_stride = survivor_stride_;
        const int m               = m_;
        const bool is_cyclic      = this->is_cyclic_;
        Vector<T> main_diagonal   = this->main_diagonal_;
        Vector<T> gamma           = gamma_;

        using MDPolicy = Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>;
        MDPolicy policy({0, 0}, {effective_batch_count, n});

        Kokkos::parallel_for(
            "SolveDiagonalCRPCR", policy, KOKKOS_LAMBDA(const int k, const int i) {
                const int batch_idx = batch_stride * k + batch_offset;
                const int offset    = batch_idx * n;

                bool is_survivor;
                if (L == 0) {
                    is_survivor = true;
                }
                else {
                    is_survivor = (i >= survivor_base) && (((i - survivor_base) % survivor_stride) == 0) &&
                                  (((i - survivor_base) / survivor_stride) < m);
                }

                const T b_i = is_survivor ? main_diagonal(offset + i) : T(1) / main_diagonal(offset + i);

                if (is_cyclic && i == 0) {
                    rhs(offset + 0) /= (b_i + gamma(batch_idx));
                }
                else {
                    rhs(offset + i) /= b_i;
                }
            });

        Kokkos::fence();
    }

private:
    static int compute_pcr_steps(int survivor_count)
    {
        int steps  = 0;
        int stride = 1;
        while (stride < survivor_count) {
            ++steps;
            // survivor_count is bounded by PCR_TARGET_SIZE, so this cannot
            // overflow for the configured target.
            stride <<= 1;
        }
        return steps;
    }

    int compute_effective_batch_count(int batch_offset, int batch_stride) const
    {
        if (batch_offset == this->batch_count_) {
            return 0;
        }
        return (this->batch_count_ - batch_offset + batch_stride - 1) / batch_stride;
    }

    static int compute_num_cr_levels(int matrix_dimension)
    {
        if (matrix_dimension <= PCR_TARGET_SIZE) {
            return 0;
        }
        int L = 1;
        while (true) {
            const long long base   = (1LL << L) - 1;
            const long long stride = 1LL << L;
            if (base >= matrix_dimension) {
                return L; // Guard; not expected to trigger for valid n.
            }
            const long long m = (matrix_dimension - 1 - base) / stride + 1;
            if (m <= PCR_TARGET_SIZE) {
                return L;
            }
            ++L;
        }
    }

    static int compute_survivor_count(int matrix_dimension, int L)
    {
        if (L == 0) {
            return matrix_dimension;
        }
        const long long base   = (1LL << L) - 1;
        const long long stride = 1LL << L;
        return static_cast<int>((matrix_dimension - 1 - base) / stride + 1);
    }

    int num_cr_levels_; // L: number of CR levels executed (0 if n <= PCR_TARGET_SIZE).
    int survivor_base_; // 2^L - 1 (0 if L==0).
    int survivor_stride_; // 2^L (1 if L==0).
    int m_; // CR survivor count / compact PCR dimension, m <= PCR_TARGET_SIZE.
    int pcr_num_steps_; // Number of compact PCR reduction steps, 0 if m<=1.

    Vector<T> q_right_trajectory_; // Persistent normalized right-coupling factors for CR-eliminated indices.

    Vector<T> cyclic_corner_; // Original cyclic corner coefficient, saved before CR overwrites sub_diagonal_.
    Vector<T> gamma_; // Sherman-Morrison correction factor, one per batch (cyclic systems only).

    Vector<T> pcr_k1_trajectory_; // Compact PCR left-elimination multipliers, O(batch*m*log m), m<=128.
    Vector<T> pcr_k2_trajectory_; // Compact PCR right-elimination multipliers, O(batch*m*log m), m<=128.

    bool is_factorized_;
};

} // namespace gmgpolar