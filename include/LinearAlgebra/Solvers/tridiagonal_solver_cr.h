#pragma once

/**
 * @brief Batched tridiagonal solver based on classical Cyclic Reduction (CR).
 *
 * Classical CR reduces the active system size by approximately a factor of
 * two per stage. Unlike Parallel Cyclic Reduction (PCR), CR does not update
 * every equation at every stage: at each level, exactly half of the
 * currently active equations are eliminated and folded into the equations
 * that survive to the next, coarser level. Reduction continues until a
 * single root equation remains.
 *
 * The solver supports arbitrary matrix dimensions without power-of-two
 * padding. Each non-root equation is eliminated exactly once, so the
 * reduction factors can be stored in O(n) additional persistent storage per
 * batch system, rather than PCR's O(n log n) trajectory.
 *
 * Each batch system is assigned one Kokkos team. The team size is selected
 * using Kokkos::AUTO; equations are distributed across team members using
 * strided ("grid-stride") loops, so the implementation does not depend on
 * the team size being equal to the matrix dimension.
 *
 * Persistent meaning of the stored arrays after setup():
 *
 *   - For every index i that is NOT the root:
 *       main_diagonal(i)      = 1 / b_i        (inverse pivot at elimination)
 *       sub_diagonal(i)       = a_i / b_i       (normalized left coupling, "qLeft")
 *       q_right_trajectory(i) = c_i / b_i       (normalized right coupling, "qRight")
 *     where a_i, b_i, c_i are the coefficients of equation i at the CR level
 *     at which i was eliminated.
 *
 *   - For the root index (root_index_):
 *       main_diagonal(root) = final reduced diagonal (NOT inverted).
 *
 * The sub-diagonal array is reused for both the "a" (left) and "c" (right)
 * coupling of the original problem via the symmetry a_i = c_{i-1}; only one
 * off-diagonal array is required, consistent with BatchedTridiagonalSolverThomas
 * and BatchedTridiagonalSolverPCR.
 *
 * For cyclic systems, the same Sherman-Morrison rank-one correction used by
 * Thomas/PCR is applied before CR factorization. Because CR overwrites the
 * sub-diagonal array in place, the original cyclic corner coefficient is
 * copied to cyclic_corner_ once, before any CR level runs, and is never
 * re-derived from sub_diagonal_ afterwards.
 *
 * setup() performs the coefficient factorization once; solve() may be
 * called any number of times afterward and never modifies the persistent
 * factorization (main_diagonal_, sub_diagonal_, q_right_trajectory_,
 * gamma_, cyclic_corner_).
 */

#include <Kokkos_Core.hpp>
#include <cstddef>
#include <stdexcept>
#include <string>

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

namespace gmgpolar
{

template <typename T>
class BatchedTridiagonalSolverCR : public BatchedTridiagonalSolverBase<T>
{
public:
    BatchedTridiagonalSolverCR(int matrix_dimension, int batch_count, bool is_cyclic = true)
        : BatchedTridiagonalSolverBase<T>(matrix_dimension, batch_count, is_cyclic)
        , num_steps_(compute_num_steps(matrix_dimension))
        , root_index_(matrix_dimension > 1 ? ((1 << num_steps_) - 1) : 0)
        , q_right_trajectory_("BatchedTridiagonalSolverCR::q_right_trajectory",
                              static_cast<std::size_t>(matrix_dimension) * static_cast<std::size_t>(batch_count))
        , cyclic_corner_("BatchedTridiagonalSolverCR::cyclic_corner", is_cyclic ? batch_count : 0)
        , gamma_("BatchedTridiagonalSolverCR::gamma", is_cyclic ? batch_count : 0)
        , is_factorized_(false)
    {
        assign(q_right_trajectory_, T(0));
        assign(cyclic_corner_, T(0));
        assign(gamma_, T(0));
    }

    /**
     * @brief Returns the stored normalized right-coupling factor q_right(i) = c_i / b_i
     *        for an equation eliminated during setup().
     *
     * Only meaningful for indices that are not the root index.
     *
     * @param batch_idx Batch index.
     * @param index     Original equation index.
     */
    KOKKOS_INLINE_FUNCTION const T& q_right(int batch_idx, int index) const
    {
        return q_right_trajectory_(static_cast<std::size_t>(batch_idx) * this->matrix_dimension_ + index);
    }

    /* -------------------------------------------------------- */
    /* Setup: Classical Cyclic Reduction coefficient elimination */
    /* -------------------------------------------------------- */
    // Performs, for cyclic systems, the Sherman-Morrison diagonal adjustment
    // (and saves the original corner coefficient), then performs the CR
    // level-by-level coefficient reduction down to a single root equation.

    void setup() override
    {
        int matrix_dimension = this->matrix_dimension_;
        int num_steps        = num_steps_;
        int root_index       = root_index_;
        bool is_cyclic       = this->is_cyclic_;

        Vector<T> main_diagonal      = this->main_diagonal_;
        Vector<T> sub_diagonal       = this->sub_diagonal_;
        Vector<T> q_right_trajectory = q_right_trajectory_;
        Vector<T> gamma              = gamma_;
        Vector<T> cyclic_corner      = cyclic_corner_;

        (void)root_index; // root requires no explicit action in setup(); it is
        // simply whichever index survives every CR level.

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;

        TeamPolicy policy(this->batch_count_, Kokkos::AUTO);

        Kokkos::parallel_for(
            "SetupCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                const int batch_idx = team_member.league_rank();
                const int offset    = batch_idx * matrix_dimension;
                const int team_size = team_member.team_size();
                const int rank      = team_member.team_rank();

                // ------------------------------------------------- //
                // Sherman-Morrison Adjustment (cyclic systems only)  //
                // - Save the original cyclic corner coefficient      //
                //   BEFORE it can be overwritten by CR.               //
                // - Modify the first and last main diagonal element. //
                // - Compute and store gamma for later use.            //
                // ------------------------------------------------- //
                if (is_cyclic) {
                    if (rank == 0) {
                        const T corner_element   = sub_diagonal(offset + matrix_dimension - 1);
                        cyclic_corner(batch_idx) = corner_element;
                        gamma(batch_idx)         = -main_diagonal(offset + 0);
                        main_diagonal(offset + 0) -= gamma(batch_idx);
                        main_diagonal(offset + matrix_dimension - 1) -=
                            corner_element * corner_element / gamma(batch_idx);
                    }
                    team_member.team_barrier();
                }

                // ------------------------------- //
                // Classical Cyclic Reduction levels //
                // ------------------------------- //
                for (int step = 0; step < num_steps; ++step) {
                    const int d = 1 << step;

                    // ---------------------------------------------------- //
                    // Phase A: eliminate equations i = (d-1) + 2*d*k.       //
                    // Read current (raw) coefficients, store normalized    //
                    // elimination factors in place.                        //
                    // ---------------------------------------------------- //
                    for (int k = rank;; k += team_size) {
                        const int i = (d - 1) + 2 * d * k;
                        if (i >= matrix_dimension) {
                            break;
                        }

                        const bool has_left  = (i - d) >= 0;
                        const bool has_right = (i + d) < matrix_dimension;

                        const T b_i = main_diagonal(offset + i);
                        const T a_i = has_left ? sub_diagonal(offset + i - d) : T(0);
                        const T c_i = has_right ? sub_diagonal(offset + i) : T(0);

                        const T inv_b  = T(1) / b_i;
                        const T qLeft  = a_i * inv_b;
                        const T qRight = c_i * inv_b;

                        main_diagonal(offset + i)      = inv_b;
                        sub_diagonal(offset + i)       = qLeft;
                        q_right_trajectory(offset + i) = qRight;
                    }

                    team_member.team_barrier();

                    // ---------------------------------------------------- //
                    // Phase B: update equations surviving to the next     //
                    // level, s = (2*d-1) + 2*d*k, using ONLY the just-     //
                    // stored normalized factors of the flanking eliminated //
                    // equations (eL = s-d always exists; eR = s+d may not).//
                    // ---------------------------------------------------- //
                    for (int k = rank;; k += team_size) {
                        const int s = (2 * d - 1) + 2 * d * k;
                        if (s >= matrix_dimension) {
                            break;
                        }

                        const int eL      = s - d;
                        const int eR      = s + d;
                        const bool has_eR = eR < matrix_dimension;

                        const T qRightEL = q_right_trajectory(offset + eL);
                        const T b_eL     = T(1) / main_diagonal(offset + eL);

                        const T c_s_raw = sub_diagonal(offset + s);
                        const T b_s_raw = main_diagonal(offset + s);

                        const T term1 = qRightEL * qRightEL * b_eL;

                        T term2 = T(0);
                        T c_new = T(0);
                        if (has_eR) {
                            const T qLeftER  = sub_diagonal(offset + eR); // already qLeft(eR) post Phase A
                            const T qRightER = q_right_trajectory(offset + eR);
                            term2            = qLeftER * c_s_raw;
                            c_new            = -c_s_raw * qRightER;
                        }

                        main_diagonal(offset + s) = b_s_raw - term1 - term2;
                        sub_diagonal(offset + s)  = c_new;
                    }

                    team_member.team_barrier();
                }
                // After the loop, main_diagonal(offset + root_index) holds the
                // final reduced (non-inverted) root diagonal.
            });

        Kokkos::fence();
        is_factorized_ = true;
    }

    /* ---------------------------------------------------- */
    /* Solve: CR forward RHS reduction + backward substitution */
    /* ---------------------------------------------------- */
    // Reuses the factorization from setup(). Eliminated-equation RHS values
    // are preserved (never overwritten) during forward reduction so they
    // remain available for backward substitution.

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

        int effective_batch_count = (this->batch_count_ - batch_offset + batch_stride - 1) / batch_stride;
        if (effective_batch_count < 0) {
            effective_batch_count = 0;
        }

        int matrix_dimension = this->matrix_dimension_;
        int num_steps        = num_steps_;
        int root_index       = root_index_;
        bool is_cyclic       = this->is_cyclic_;

        Vector<T> main_diagonal      = this->main_diagonal_;
        Vector<T> sub_diagonal       = this->sub_diagonal_;
        Vector<T> q_right_trajectory = q_right_trajectory_;
        Vector<T> gamma              = gamma_;
        Vector<T> cyclic_corner      = cyclic_corner_;

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;

        if (!is_cyclic) {
            TeamPolicy policy(effective_batch_count, Kokkos::AUTO);

            Kokkos::parallel_for(
                "SolveCRNonCyclic", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                    const int k         = team_member.league_rank();
                    const int batch_idx = batch_stride * k + batch_offset;
                    const int offset    = batch_idx * matrix_dimension;
                    const int team_size = team_member.team_size();
                    const int rank      = team_member.team_rank();

                    // -------------------- //
                    // Forward RHS Reduction //
                    // -------------------- //
                    for (int step = 0; step < num_steps; ++step) {
                        const int d = 1 << step;

                        for (int kk = rank;; kk += team_size) {
                            const int s = (2 * d - 1) + 2 * d * kk;
                            if (s >= matrix_dimension) {
                                break;
                            }

                            const int eL      = s - d;
                            const int eR      = s + d;
                            const bool has_eR = eR < matrix_dimension;

                            const T factor_left  = q_right_trajectory(offset + eL);
                            const T factor_right = has_eR ? sub_diagonal(offset + eR) : T(0);
                            const T rhs_eR       = has_eR ? rhs(offset + eR) : T(0);

                            rhs(offset + s) -= factor_left * rhs(offset + eL) + factor_right * rhs_eR;
                        }

                        team_member.team_barrier();
                    }

                    // ----------- //
                    // Root Solve  //
                    // ----------- //
                    if (rank == 0) {
                        rhs(offset + root_index) /= main_diagonal(offset + root_index);
                    }
                    team_member.team_barrier();

                    // ---------------------- //
                    // Backward Substitution  //
                    // ---------------------- //
                    for (int step = num_steps - 1; step >= 0; --step) {
                        const int d = 1 << step;

                        for (int kk = rank;; kk += team_size) {
                            const int i = (d - 1) + 2 * d * kk;
                            if (i >= matrix_dimension) {
                                break;
                            }

                            const bool has_left  = (i - d) >= 0;
                            const bool has_right = (i + d) < matrix_dimension;

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
            // Cyclic solve: reduce the RHS and the Sherman-Morrison auxiliary
            // vector simultaneously, using the same stored CR factors. The
            // auxiliary vector uses per-team scratch memory rather than a
            // second permanent batch_count * n allocation.
            const std::size_t scratch_bytes = static_cast<std::size_t>(matrix_dimension) * sizeof(T);

            TeamPolicy policy(effective_batch_count, Kokkos::AUTO);
            policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

            Kokkos::parallel_for(
                "SolveCRCyclic", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                    const int k         = team_member.league_rank();
                    const int batch_idx = batch_stride * k + batch_offset;
                    const int offset    = batch_idx * matrix_dimension;
                    const int team_size = team_member.team_size();
                    const int rank      = team_member.team_rank();

                    T* buffer = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));

                    // ------------------------------------------------ //
                    // Initialize the Sherman-Morrison auxiliary vector: //
                    // u = gamma * e_0 + corner * e_{n-1}                //
                    // ------------------------------------------------ //
                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        buffer[i] = T(0);
                    }
                    team_member.team_barrier();
                    if (rank == 0) {
                        buffer[0] = gamma(batch_idx);
                        if (matrix_dimension > 1) {
                            buffer[matrix_dimension - 1] = cyclic_corner(batch_idx);
                        }
                    }
                    team_member.team_barrier();

                    // -------------------- //
                    // Forward RHS Reduction //
                    // -------------------- //
                    for (int step = 0; step < num_steps; ++step) {
                        const int d = 1 << step;

                        for (int kk = rank;; kk += team_size) {
                            const int s = (2 * d - 1) + 2 * d * kk;
                            if (s >= matrix_dimension) {
                                break;
                            }

                            const int eL      = s - d;
                            const int eR      = s + d;
                            const bool has_eR = eR < matrix_dimension;

                            const T factor_left  = q_right_trajectory(offset + eL);
                            const T factor_right = has_eR ? sub_diagonal(offset + eR) : T(0);
                            const T rhs_eR       = has_eR ? rhs(offset + eR) : T(0);
                            const T buf_eR       = has_eR ? buffer[eR] : T(0);

                            rhs(offset + s) -= factor_left * rhs(offset + eL) + factor_right * rhs_eR;
                            buffer[s] -= factor_left * buffer[eL] + factor_right * buf_eR;
                        }

                        team_member.team_barrier();
                    }

                    // ----------- //
                    // Root Solve  //
                    // ----------- //
                    if (rank == 0) {
                        rhs(offset + root_index) /= main_diagonal(offset + root_index);
                        buffer[root_index] /= main_diagonal(offset + root_index);
                    }
                    team_member.team_barrier();

                    // ---------------------- //
                    // Backward Substitution  //
                    // ---------------------- //
                    for (int step = num_steps - 1; step >= 0; --step) {
                        const int d = 1 << step;

                        for (int kk = rank;; kk += team_size) {
                            const int i = (d - 1) + 2 * d * kk;
                            if (i >= matrix_dimension) {
                                break;
                            }

                            const bool has_left  = (i - d) >= 0;
                            const bool has_right = (i + d) < matrix_dimension;

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

                    // ------------------------------- //
                    // Sherman-Morrison Reconstruction //
                    // ------------------------------- //
                    const T corner = cyclic_corner(batch_idx);
                    const T g      = gamma(batch_idx);

                    const T dot_product_x_v = rhs(offset + 0) + corner / g * rhs(offset + matrix_dimension - 1);
                    const T dot_product_u_v = buffer[0] + corner / g * buffer[matrix_dimension - 1];
                    const T factor          = dot_product_x_v / (T(1) + dot_product_u_v);

                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        rhs(offset + i) -= factor * buffer[i];
                    }
                });
        }
        Kokkos::fence();
    }

    /* ---------------------------- */
    /* Solve: Diagonal Scaling Only */
    /* ---------------------------- */
    // Valid when the underlying matrix has zero off-diagonal coupling. Does
    // NOT perform CR reduction. Because CR's persistent storage inverts the
    // diagonal at every eliminated (non-root) index, the true diagonal b_i is
    // reconstructed as 1/main_diagonal(i) there, and used directly (raw) at
    // the root index. For cyclic systems, the same gamma correction used by
    // Thomas/PCR is applied to index 0.

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

        int effective_batch_count = (this->batch_count_ - batch_offset + batch_stride - 1) / batch_stride;
        if (effective_batch_count < 0) {
            effective_batch_count = 0;
        }

        const int matrix_dimension = this->matrix_dimension_;
        const int root_index       = root_index_;
        const bool is_cyclic       = this->is_cyclic_;
        Vector<T> main_diagonal    = this->main_diagonal_;
        Vector<T> gamma            = gamma_;

        using MDPolicy = Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>;
        MDPolicy policy({0, 0}, {effective_batch_count, matrix_dimension});

        Kokkos::parallel_for(
            "SolveDiagonalCR", policy, KOKKOS_LAMBDA(const int k, const int i) {
                const int batch_idx = batch_stride * k + batch_offset;
                const int offset    = batch_idx * matrix_dimension;

                const T b_i = (i == root_index) ? main_diagonal(offset + i) : T(1) / main_diagonal(offset + i);

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
    static int compute_num_steps(int matrix_dimension)
    {
        // floor(log2(matrix_dimension)) for matrix_dimension > 1, else 0.
        // Uses integer bit logic only (no std::log2), per repository policy.
        if (matrix_dimension <= 1) {
            return 0;
        }
        int steps = 0;
        while ((1 << (steps + 1)) <= matrix_dimension) {
            ++steps;
        }
        return steps;
    }

    int num_steps_; // Number of CR reduction levels: floor(log2(n)), 0 for n<=1.
    int root_index_; // Index of the single equation surviving all CR levels.

    Vector<T> q_right_trajectory_; // Persistent normalized right-coupling factors, O(batch_count * n).

    Vector<T> cyclic_corner_; // Original cyclic corner coefficient, saved before CR overwrites sub_diagonal_.
    Vector<T> gamma_; // Sherman-Morrison correction factor, one per batch (cyclic systems only).

    bool is_factorized_;
};

} // namespace gmgpolar