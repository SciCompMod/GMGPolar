#pragma once

/**
 * @brief Batched tridiagonal solver based on Parallel Cyclic Reduction (PCR).
 *
 * The solver operates on a batch of independent tridiagonal systems and uses
 * Kokkos teams to distribute the equations of each system across team members.
 *
 * PCR performs O(n log n) arithmetic work with O(log n) parallel depth. This
 * allows the equations within a system to be processed concurrently, which is
 * beneficial when the number of systems is small relative to the system size.
 *
 * Each system is assigned one Kokkos team. The team size is selected using
 * Kokkos::AUTO to allow the execution backend to choose an appropriate value.
 * The equations are distributed across team members using strided loops, so
 * the implementation does not depend on the team size being equal to the
 * matrix dimension.
 *
 * @tparam T Scalar type used for matrix coefficients and right-hand sides.
 *
 * @note setup() performs the coefficient reduction and stores the reduction
 *       coefficients. Subsequent calls to solve() reuse this factorization.
 *
 * For cyclic systems, setup() also prepares the diagonal correction required
 * by the Sherman–Morrison reconstruction. The coefficient reduction is shared
 * by the right-hand side and auxiliary solve performed by solve().
 */

#include <Kokkos_Core.hpp>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

namespace gmgpolar
{

/**
 * @brief Computes the clamped left and right neighbor indices used by PCR.
 *
 * The boundary convention is that the left coefficient of the first equation
 * and the right coefficient of the last equation are zero. Clamping the
 * neighbor indices to the valid range allows the same PCR update to be used
 * for boundary and interior equations.
 *
 * @param i Current equation index.
 * @param delta Distance to the neighboring equation for the current PCR step.
 * @param n Number of equations.
 * @param[out] iLeft Left neighbor index.
 * @param[out] iRight Right neighbor index.
 */
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
class BatchedTridiagonalSolverPCR : public BatchedTridiagonalSolverBase<T>
{
public:
    BatchedTridiagonalSolverPCR(int matrix_dimension, int batch_count, bool is_cyclic = true)
        : BatchedTridiagonalSolverBase<T>(matrix_dimension, batch_count, is_cyclic)
        , gamma_("BatchedTridiagonalSolverPCR::gamma", is_cyclic ? batch_count : 0)
        , is_factorized_(false)
        , num_steps_(
              matrix_dimension > 1 ? static_cast<int>(std::ceil(std::log2(static_cast<double>(matrix_dimension)))) : 0)
        , k1_trajectory_("BatchedTridiagonalSolverPCR::k1_trajectory", static_cast<std::size_t>(batch_count) *
                                                                           static_cast<std::size_t>(num_steps_) *
                                                                           static_cast<std::size_t>(matrix_dimension))
        , k2_trajectory_("BatchedTridiagonalSolverPCR::k2_trajectory", static_cast<std::size_t>(batch_count) *
                                                                           static_cast<std::size_t>(num_steps_) *
                                                                           static_cast<std::size_t>(matrix_dimension))
    {
        assign(gamma_, T(0));
        assign(k1_trajectory_, T(0));
        assign(k2_trajectory_, T(0));
    }

    /**
     * @brief Returns a stored PCR left reduction coefficient.
     *
     * @param batch_idx Batch index.
     * @param step PCR reduction step.
     * @param index Equation index.
     */
    KOKKOS_INLINE_FUNCTION const T& k1(const int batch_idx, const int step, const int index) const
    {
        return k1_trajectory_(static_cast<std::size_t>(batch_idx) * num_steps_ * this->matrix_dimension_ +
                              static_cast<std::size_t>(step) * this->matrix_dimension_ + index);
    }

    /**
     * @brief Returns a stored PCR right reduction coefficient.
     *
     * @param batch_idx Batch index.
     * @param step PCR reduction step.
     * @param index Equation index.
     */
    KOKKOS_INLINE_FUNCTION const T& k2(const int batch_idx, const int step, const int index) const
    {
        return k2_trajectory_(static_cast<std::size_t>(batch_idx) * num_steps_ * this->matrix_dimension_ +
                              static_cast<std::size_t>(step) * this->matrix_dimension_ + index);
    }

    /**
     * @brief Performs PCR coefficient reduction and prepares the solver.
     *
     * The coefficient reduction is performed once and the resulting reduction
     * coefficients are stored for reuse by solve(). The final reduced diagonal
     * replaces the original main diagonal.
     *
     * For cyclic systems, the diagonal is modified as part of the
     * Sherman–Morrison formulation and the corresponding correction factor is
     * stored in gamma_.
     *
     * The sub-diagonal remains unchanged by this operation.
     */
    void setup() override
    {
        int matrix_dimension = this->matrix_dimension_;
        int num_steps        = num_steps_;
        bool is_cyclic       = this->is_cyclic_;

        Vector<T> main_diagonal = this->main_diagonal_;
        Vector<T> sub_diagonal  = this->sub_diagonal_;
        Vector<T> gamma         = gamma_;
        Vector<T> k1_trajectory = k1_trajectory_;
        Vector<T> k2_trajectory = k2_trajectory_;

        if (matrix_dimension == 1) {
            is_factorized_ = true;
            return;
        }

        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using TeamMember = typename TeamPolicy::member_type;

        // For a symmetric tridiagonal system, the left coefficient of an
        // equation is the right coefficient of the corresponding left
        // neighbor. The left coefficients can therefore be reconstructed from
        // e[] rather than stored separately, reducing team scratch storage.
        const std::size_t scratch_bytes = 2ull * 2ull * static_cast<std::size_t>(matrix_dimension) * sizeof(T);

        TeamPolicy policy(this->batch_count_, Kokkos::AUTO);
        policy.set_scratch_size(0, Kokkos::PerTeam(static_cast<int>(scratch_bytes)));

        Kokkos::parallel_for(
            "SetupPCR", policy, KOKKOS_LAMBDA(const TeamMember& team_member) {
                const int batch_idx = team_member.league_rank();
                const int offset    = batch_idx * matrix_dimension;
                const int team_size = team_member.team_size();
                const int rank      = team_member.team_rank();

                T* scratch = static_cast<T*>(team_member.team_scratch(0).get_shmem(scratch_bytes));
                // Layout: [e0 | b0 | e1 | b1].
                T* e[2] = {scratch, scratch + 2 * matrix_dimension};
                T* b[2] = {scratch + matrix_dimension, scratch + 3 * matrix_dimension};

                int cur = 0;

                // e[] stores the right coefficients of the current reduced
                // system. The left coefficients are reconstructed from e[].
                for (int i = rank; i < matrix_dimension; i += team_size) {
                    e[cur][i] = (i == matrix_dimension - 1) ? T(0) : sub_diagonal(offset + i);
                    b[cur][i] = main_diagonal(offset + i);
                }

                team_member.team_barrier();

                if (is_cyclic) {
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

                        const T a_i      = (i >= delta) ? e[cur][i - delta] : T(0);
                        const T a_iRight = (iRight >= delta) ? e[cur][iRight - delta] : T(0);
                        const T c_i      = e[cur][i];
                        const T c_iLeft  = e[cur][iLeft];

                        const T k1_val = a_i / b[cur][iLeft];
                        const T k2_val = c_i / b[cur][iRight];

                        k1_trajectory(static_cast<std::size_t>(batch_idx) * num_steps * matrix_dimension +
                                      static_cast<std::size_t>(step) * matrix_dimension + i) = k1_val;
                        k2_trajectory(static_cast<std::size_t>(batch_idx) * num_steps * matrix_dimension +
                                      static_cast<std::size_t>(step) * matrix_dimension + i) = k2_val;

                        const int nxt = 1 - cur;
                        e[nxt][i]     = -e[cur][iRight] * k2_val;
                        b[nxt][i]     = b[cur][i] - c_iLeft * k1_val - a_iRight * k2_val;
                    }

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

    /**
     * @brief Solves the factored tridiagonal systems for the supplied RHS.
     *
     * The coefficient reduction performed by setup() is reused. The stored PCR
     * trajectory is applied to the right-hand side, followed by division by the
     * reduced diagonal.
     *
     * For cyclic systems, the right-hand side and the auxiliary vector used by
     * the Sherman–Morrison reconstruction are reduced using the same stored
     * trajectory within a single kernel launch.
     *
     * @param rhs Right-hand sides to solve. The vector is overwritten with the
     *            corresponding solutions.
     * @param batch_offset First batch index to process.
     * @param batch_stride Distance between processed batch indices.
     *
     * @throws std::runtime_error if setup() has not been called.
     */
    void solve(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1) override
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        const int effective_batch_count = (this->batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        int matrix_dimension = this->matrix_dimension_;
        int num_steps        = num_steps_;
        bool is_cyclic       = this->is_cyclic_;

        Vector<T> main_diagonal = this->main_diagonal_;
        Vector<T> sub_diagonal  = this->sub_diagonal_;
        Vector<T> gamma         = gamma_;
        Vector<T> k1_trajectory = k1_trajectory_;
        Vector<T> k2_trajectory = k2_trajectory_;

        if (matrix_dimension == 1) {
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

                        team_member.team_barrier();
                        cur = 1 - cur;
                    }

                    for (int i = rank; i < matrix_dimension; i += team_size) {
                        rhs(offset + i) = d[cur][i] / main_diagonal(offset + i);
                    }
                });
        }
        else {
            // The cyclic solve simultaneously reduces the right-hand side and
            // the auxiliary Sherman–Morrison vector. The original corner
            // coefficient remains available in this->sub_diagonal_.
            const std::size_t scratch_bytes = 4ull * static_cast<std::size_t>(matrix_dimension) * sizeof(T);

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
                    // Layout: [d_rhs(0) | d_buf(0) | d_rhs(1) | d_buf(1)].
                    T* d_rhs[2] = {scratch, scratch + 2 * matrix_dimension};
                    T* d_buf[2] = {scratch + matrix_dimension, scratch + 3 * matrix_dimension};

                    const T cyclic_corner_element = sub_diagonal(offset + matrix_dimension - 1);

                    int cur = 0;
                    // Initial auxiliary vector for the Sherman–Morrison
                    // reconstruction.
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

                        team_member.team_barrier();
                        cur = 1 - cur;
                    }

                    // The unused buffer stores the reduced solutions so all
                    // team members can access the entries required for the
                    // Sherman–Morrison reconstruction.
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

    /**
     * @brief Solves systems whose matrix has already been reduced to diagonal form.
     *
     * Each matrix entry is independent, so the operation is parallelized over
     * both batch and equation indices.
     *
     * For cyclic systems, the first diagonal entry includes the corresponding
     * Sherman–Morrison diagonal correction stored in gamma_.
     *
     * @param rhs Right-hand sides to solve. The vector is overwritten with the
     *            resulting solution.
     * @param batch_offset First batch index to process.
     * @param batch_stride Distance between processed batch indices.
     *
     * @throws std::runtime_error if setup() has not been called.
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
            "SolveDiagonal", policy, KOKKOS_LAMBDA(const int k, const int i) {
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

private:
    int num_steps_; // Number of PCR reduction steps.
    Vector<T> k1_trajectory_; // Stored left reduction coefficients.
    Vector<T> k2_trajectory_; // Stored right reduction coefficients.

    Vector<T> gamma_;
    bool is_factorized_;
};

} // namespace gmgpolar
