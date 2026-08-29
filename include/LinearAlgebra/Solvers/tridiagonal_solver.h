#pragma once

#include <Kokkos_Core.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

namespace gmgpolar
{

namespace detail
{

KOKKOS_INLINE_FUNCTION
void pcr_neighbors(
const int i,
const int delta,
const int n,
int& iLeft,
int& iRight)
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

} // namespace detail

template <typename T>
class BatchedTridiagonalSolver
{
private:
using execution_space = Kokkos::DefaultExecutionSpace;
using team_policy = Kokkos::TeamPolicy<execution_space>;
using team_member = typename team_policy::member_type;
using scratch_space = typename team_member::scratch_memory_space;


using unmanaged_scratch_view =
    Kokkos::View<
        T*,
        scratch_space,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

struct TeamSizeFunctor
{
    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member&) const
    {
    }
};

static int compute_num_steps(const int n)
{
    if (n <= 1) {
        return 0;
    }

    return static_cast<int>(
        std::ceil(std::log2(static_cast<double>(n))));
}

static std::size_t coefficient_scratch_bytes(const int n)
{
    return static_cast<std::size_t>(6) *
           static_cast<std::size_t>(n) *
           sizeof(T);
}

static std::size_t solve_scratch_bytes(
    const int n,
    const bool cyclic)
{
    const int number_of_arrays = cyclic ? 4 : 2;

    return static_cast<std::size_t>(number_of_arrays) *
           static_cast<std::size_t>(n) *
           sizeof(T);
}


public:
BatchedTridiagonalSolver(
int matrix_dimension,
int batch_count,
bool is_cyclic = true)
: matrix_dimension_(matrix_dimension)
, batch_count_(batch_count)
, num_steps_(compute_num_steps(matrix_dimension))
, main_diagonal_(
"BatchedTridiagonalSolver::main_diagonal",
matrix_dimension * batch_count)
, sub_diagonal_(
"BatchedTridiagonalSolver::sub_diagonal",
matrix_dimension * batch_count)
, buffer_(
"BatchedTridiagonalSolver::buffer",
is_cyclic ? matrix_dimension * batch_count : 0)
, gamma_(
"BatchedTridiagonalSolver::gamma",
is_cyclic ? batch_count : 0)
, k1_trajectory_(
"BatchedTridiagonalSolver::k1_trajectory",
static_cast[std::size_t](std::size_t)(batch_count) *
static_cast[std::size_t](std::size_t)(num_steps_) *
static_cast[std::size_t](std::size_t)(matrix_dimension))
, k2_trajectory_(
"BatchedTridiagonalSolver::k2_trajectory",
static_cast[std::size_t](std::size_t)(batch_count) *
static_cast[std::size_t](std::size_t)(num_steps_) *
static_cast[std::size_t](std::size_t)(matrix_dimension))
, is_cyclic_(is_cyclic)
, is_factorized_(false)
{
if (matrix_dimension_ <= 0) {
throw std::invalid_argument(
"BatchedTridiagonalSolver: "
"matrix_dimension must be positive.");
}


    if (batch_count_ < 0) {
        throw std::invalid_argument(
            "BatchedTridiagonalSolver: "
            "batch_count must be non-negative.");
    }

    /*
     * PCR uses one team thread per equation.
     *
     * Explicitly reject matrices larger than the backend-supported
     * maximum team size rather than silently falling back to a
     * slower one-thread-per-system implementation.
     */
    team_policy probe_policy(1, Kokkos::AUTO);
    const int max_team_size =
        probe_policy.team_size_max(TeamSizeFunctor());

    if (matrix_dimension_ > max_team_size) {
        throw std::runtime_error(
            "BatchedTridiagonalSolver: matrix dimension " +
            std::to_string(matrix_dimension_) +
            " exceeds the Kokkos backend maximum team size " +
            std::to_string(max_team_size) +
            ". PCR requires one team thread per equation.");
    }

    assign(main_diagonal_, T(0));
    assign(sub_diagonal_, T(0));
    assign(k1_trajectory_, T(0));
    assign(k2_trajectory_, T(0));
}

/* ------------------- */
/* Accessors for sizes */
/* ------------------- */

KOKKOS_INLINE_FUNCTION
int matrixDimension() const
{
    return matrix_dimension_;
}

KOKKOS_INLINE_FUNCTION
int batchCount() const
{
    return batch_count_;
}

/* ---------------------------- */
/* Accessors for matrix entries */
/* ---------------------------- */

KOKKOS_INLINE_FUNCTION
const T& main_diagonal(
    const int batch_idx,
    const int index) const
{
    return main_diagonal_(
        batch_idx * matrix_dimension_ + index);
}

KOKKOS_INLINE_FUNCTION
void set_main_diagonal(
    const int batch_idx,
    const int index,
    const T& value) const
{
    main_diagonal_(
        batch_idx * matrix_dimension_ + index) = value;
}

KOKKOS_INLINE_FUNCTION
void increase_main_diagonal(
    const int batch_idx,
    const int index,
    const T& value) const
{
    main_diagonal_(
        batch_idx * matrix_dimension_ + index) += value;
}

KOKKOS_INLINE_FUNCTION
const T& sub_diagonal(
    const int batch_idx,
    const int index) const
{
    return sub_diagonal_(
        batch_idx * matrix_dimension_ + index);
}

KOKKOS_INLINE_FUNCTION
void set_sub_diagonal(
    const int batch_idx,
    const int index,
    const T& value) const
{
    sub_diagonal_(
        batch_idx * matrix_dimension_ + index) = value;
}

KOKKOS_INLINE_FUNCTION
void increase_sub_diagonal(
    const int batch_idx,
    const int index,
    const T& value) const
{
    sub_diagonal_(
        batch_idx * matrix_dimension_ + index) += value;
}

KOKKOS_INLINE_FUNCTION
const T& cyclic_corner(const int batch_idx) const
{
    return sub_diagonal_(
        batch_idx * matrix_dimension_ +
        (matrix_dimension_ - 1));
}

KOKKOS_INLINE_FUNCTION
T& set_cyclic_corner(
    const int batch_idx,
    const T& value) const
{
    return sub_diagonal_(
        batch_idx * matrix_dimension_ +
        (matrix_dimension_ - 1)) = value;
}

KOKKOS_INLINE_FUNCTION
void increase_cyclic_corner(
    const int batch_idx,
    const T& value) const
{
    sub_diagonal_(
        batch_idx * matrix_dimension_ +
        (matrix_dimension_ - 1)) += value;
}

/* ---------------------------------------------- */
/* Setup: Parallel Cyclic Reduction                */
/* ---------------------------------------------- */

void setup()
{
    const int n = matrix_dimension_;
    const int num_steps = num_steps_;
    const int batch_count = batch_count_;
    const bool cyclic = is_cyclic_;

    Vector<T> main_diagonal = main_diagonal_;
    Vector<T> sub_diagonal = sub_diagonal_;
    Vector<T> gamma = gamma_;
    Vector<T> k1_trajectory = k1_trajectory_;
    Vector<T> k2_trajectory = k2_trajectory_;

    /*
     * Six arrays are required:
     *
     *   a0, b0, c0
     *   a1, b1, c1
     *
     * The two sets are ping-pong buffers for PCR.
     */
    const std::size_t scratch_bytes =
        coefficient_scratch_bytes(n);

    team_policy policy(batch_count, n);

    policy.set_scratch_size(
        0,
        Kokkos::PerTeam(scratch_bytes));

    Kokkos::parallel_for(
        "BatchedTridiagonalSolver::PCRSetup",
        policy,
        KOKKOS_LAMBDA(const team_member& team)
        {
            const int batch_idx = team.league_rank();
            const int i = team.team_rank();
            const int offset = batch_idx * n;

            /*
             * Allocate each scratch array consecutively.
             *
             * Calling get_shmem() advances the team's scratch
             * allocator. This is important: constructing every
             * unmanaged View directly from team_scratch(0) would
             * make all of the Views alias the same memory.
             */
            unmanaged_scratch_view a0(
                team.team_scratch(0).get_shmem(
                    static_cast<std::size_t>(n) * sizeof(T)),
                n);

            unmanaged_scratch_view b0(
                team.team_scratch(0).get_shmem(
                    static_cast<std::size_t>(n) * sizeof(T)),
                n);

            unmanaged_scratch_view c0(
                team.team_scratch(0).get_shmem(
                    static_cast<std::size_t>(n) * sizeof(T)),
                n);

            unmanaged_scratch_view a1(
                team.team_scratch(0).get_shmem(
                    static_cast<std::size_t>(n) * sizeof(T)),
                n);

            unmanaged_scratch_view b1(
                team.team_scratch(0).get_shmem(
                    static_cast<std::size_t>(n) * sizeof(T)),
                n);

            unmanaged_scratch_view c1(
                team.team_scratch(0).get_shmem(
                    static_cast<std::size_t>(n) * sizeof(T)),
                n);

            /*
             * Initialize the original tridiagonal system.
             *
             * The boundary entries a[0] and c[n-1] are virtual zero
             * coefficients.
             */
            a0(i) =
                (i == 0)
                    ? T(0)
                    : sub_diagonal(offset + i - 1);

            b0(i) =
                main_diagonal(offset + i);

            c0(i) =
                (i == n - 1)
                    ? T(0)
                    : sub_diagonal(offset + i);

            /*
             * Preserve sub_diagonal_ exactly as supplied by the caller.
             *
             * The old Thomas implementation overwrote sub_diagonal_
             * with elimination multipliers. PCR does not do that:
             * its multipliers are stored in k1_trajectory_ and
             * k2_trajectory_. Keeping the original coefficients makes
             * it possible to call setup() again after matrix changes.
             *
             * For n == 1 there is no genuine cyclic off-diagonal
             * connection, so the scalar system is handled as a
             * one-equation system without a Sherman-Morrison update.
             */
            if (cyclic && n > 1 && i == 0) {
                const T cyclic_corner_element =
                    sub_diagonal(offset + n - 1);

                gamma(batch_idx) =
                    -main_diagonal(offset);

                b0(0) -= gamma(batch_idx);

                b0(n - 1) -=
                    cyclic_corner_element *
                    cyclic_corner_element /
                    gamma(batch_idx);
            }

            team.team_barrier();

            int current = 0;

            for (int step = 0;
                 step < num_steps;
                 ++step) {
                const int delta = 1 << step;

                unmanaged_scratch_view a_current =
                    (current == 0) ? a0 : a1;

                unmanaged_scratch_view b_current =
                    (current == 0) ? b0 : b1;

                unmanaged_scratch_view c_current =
                    (current == 0) ? c0 : c1;

                unmanaged_scratch_view a_next =
                    (current == 0) ? a1 : a0;

                unmanaged_scratch_view b_next =
                    (current == 0) ? b1 : b0;

                unmanaged_scratch_view c_next =
                    (current == 0) ? c1 : c0;

                int iLeft;
                int iRight;

                detail::pcr_neighbors(
                    i,
                    delta,
                    n,
                    iLeft,
                    iRight);

                const T k1 =
                    a_current(i) /
                    b_current(iLeft);

                const T k2 =
                    c_current(i) /
                    b_current(iRight);

                const std::size_t trajectory_index =
                    static_cast<std::size_t>(batch_idx) *
                        static_cast<std::size_t>(num_steps) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(step) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(i);

                k1_trajectory(trajectory_index) = k1;
                k2_trajectory(trajectory_index) = k2;

                a_next(i) =
                    -a_current(iLeft) * k1;

                b_next(i) =
                    b_current(i) -
                    c_current(iLeft) * k1 -
                    a_current(iRight) * k2;

                c_next(i) =
                    -c_current(iRight) * k2;

                /*
                 * PCR reads neighbors' pre-update values.
                 * Do not swap the ping-pong buffers until every
                 * team thread has completed this step.
                 */
                team.team_barrier();

                current = 1 - current;
            }

            unmanaged_scratch_view b_final =
                (current == 0) ? b0 : b1;

            main_diagonal(offset + i) =
                b_final(i);
        });

    Kokkos::fence();

    is_factorized_ = true;
}

/* ---------------------------------------- */
/* Solve: Parallel Cyclic Reduction         */
/* ---------------------------------------- */

void solve(
    Vector<T> rhs,
    int batch_offset = 0,
    int batch_stride = 1)
{
    if (!is_factorized_) {
        throw std::runtime_error(
            "Error: Matrix must be factorized before solving.");
    }

    const int effective_batch_count =
        (batch_count_ - batch_offset + batch_stride - 1) /
        batch_stride;

    const int n = matrix_dimension_;
    const int num_steps = num_steps_;
    const int offset_batch = batch_offset;
    const int stride = batch_stride;
    const bool cyclic = is_cyclic_;

    Vector<T> main_diagonal = main_diagonal_;
    Vector<T> sub_diagonal = sub_diagonal_;
    Vector<T> gamma = gamma_;
    Vector<T> k1_trajectory = k1_trajectory_;
    Vector<T> k2_trajectory = k2_trajectory_;

    /*
     * Non-cyclic:
     *
     *   d_rhs0, d_rhs1
     *
     * Cyclic:
     *
     *   d_rhs0, d_rhs1
     *   d_buf0, d_buf1
     */
    const std::size_t scratch_bytes =
        solve_scratch_bytes(n, cyclic);

    team_policy policy(
        effective_batch_count,
        n);

    policy.set_scratch_size(
        0,
        Kokkos::PerTeam(scratch_bytes));

    Kokkos::parallel_for(
        "BatchedTridiagonalSolver::PCRSolve",
        policy,
        KOKKOS_LAMBDA(const team_member& team)
        {
            const int k = team.league_rank();
            const int batch_idx =
                stride * k + offset_batch;
            const int offset = batch_idx * n;
            const int i = team.team_rank();

            /*
             * Allocate RHS ping-pong buffers.
             */
            unmanaged_scratch_view d_rhs0(
                team.team_scratch(0).get_shmem(
                    static_cast<std::size_t>(n) * sizeof(T)),
                n);

            unmanaged_scratch_view d_rhs1(
                team.team_scratch(0).get_shmem(
                    static_cast<std::size_t>(n) * sizeof(T)),
                n);

            /*
             * Allocate cyclic auxiliary-RHS ping-pong buffers.
             */
            unmanaged_scratch_view d_buf0;
            unmanaged_scratch_view d_buf1;

            if (cyclic) {
                d_buf0 = unmanaged_scratch_view(
                    team.team_scratch(0).get_shmem(
                        static_cast<std::size_t>(n) * sizeof(T)),
                    n);

                d_buf1 = unmanaged_scratch_view(
                    team.team_scratch(0).get_shmem(
                        static_cast<std::size_t>(n) * sizeof(T)),
                    n);
            }

            d_rhs0(i) =
                rhs(offset + i);

            /*
             * Initial Sherman-Morrison auxiliary RHS.
             *
             * The PCR reduction itself performs all elimination
             * transformations, so unlike the Thomas implementation
             * there is no sequential chain multiplication here.
             */
            if (cyclic) {
                const T cyclic_corner_element =
                    sub_diagonal(offset + n - 1);

                if (n > 1) {
                    d_buf0(i) =
                        (i == 0)
                            ? gamma(batch_idx)
                            : ((i == n - 1)
                                   ? cyclic_corner_element
                                   : T(0));
                }
                else {
                    d_buf0(i) = T(0);
                }
            }

            team.team_barrier();

            int current = 0;

            for (int step = 0;
                 step < num_steps;
                 ++step) {
                const int delta = 1 << step;

                unmanaged_scratch_view d_rhs_current =
                    (current == 0)
                        ? d_rhs0
                        : d_rhs1;

                unmanaged_scratch_view d_rhs_next =
                    (current == 0)
                        ? d_rhs1
                        : d_rhs0;

                unmanaged_scratch_view d_buf_current;
                unmanaged_scratch_view d_buf_next;

                if (cyclic) {
                    d_buf_current =
                        (current == 0)
                            ? d_buf0
                            : d_buf1;

                    d_buf_next =
                        (current == 0)
                            ? d_buf1
                            : d_buf0;
                }

                int iLeft;
                int iRight;

                detail::pcr_neighbors(
                    i,
                    delta,
                    n,
                    iLeft,
                    iRight);

                const std::size_t trajectory_index =
                    static_cast<std::size_t>(batch_idx) *
                        static_cast<std::size_t>(num_steps) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(step) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(i);

                const T k1 =
                    k1_trajectory(trajectory_index);

                const T k2 =
                    k2_trajectory(trajectory_index);

                d_rhs_next(i) =
                    d_rhs_current(i) -
                    d_rhs_current(iLeft) * k1 -
                    d_rhs_current(iRight) * k2;

                if (cyclic) {
                    d_buf_next(i) =
                        d_buf_current(i) -
                        d_buf_current(iLeft) * k1 -
                        d_buf_current(iRight) * k2;
                }

                /*
                 * PCR reads neighbors' pre-update values.
                 * The barrier must occur before swapping the
                 * ping-pong buffers.
                 */
                team.team_barrier();

                current = 1 - current;
            }

            unmanaged_scratch_view d_rhs_final =
                (current == 0)
                    ? d_rhs0
                    : d_rhs1;

            d_rhs_final(i) /=
                main_diagonal(offset + i);

            if (!cyclic || n == 1) {
                rhs(offset + i) =
                    d_rhs_final(i);

                return;
            }

            unmanaged_scratch_view d_buf_final =
                (current == 0)
                    ? d_buf0
                    : d_buf1;

            d_buf_final(i) /=
                main_diagonal(offset + i);

            team.team_barrier();

            const T cyclic_corner_element =
                sub_diagonal(offset + n - 1);

            T factor = T(0);

            if (team.team_rank() == 0) {
                const T dot_product_x_v =
                    d_rhs_final(0) +
                    (cyclic_corner_element /
                     gamma(batch_idx)) *
                        d_rhs_final(n - 1);

                const T dot_product_u_v =
                    d_buf_final(0) +
                    (cyclic_corner_element /
                     gamma(batch_idx)) *
                        d_buf_final(n - 1);

                factor =
                    dot_product_x_v /
                    (T(1) + dot_product_u_v);
            }

            team.team_broadcast(
                factor,
                0);

            rhs(offset + i) =
                d_rhs_final(i) -
                factor * d_buf_final(i);
        });

    Kokkos::fence();
}

/* ---------------------------- */
/* Solve: Diagonal Scaling Only */
/* ---------------------------- */

void solve_diagonal(
    Vector<T> rhs,
    int batch_offset = 0,
    int batch_stride = 1)
{
    if (!is_factorized_) {
        throw std::runtime_error(
            "Error: Matrix must be factorized before solving.");
    }

    // Compute the effective number of batches to solve
    int effective_batch_count =
        (batch_count_ - batch_offset + batch_stride - 1) /
        batch_stride;

    // Create local copies for lambda capture
    int matrix_dimension = matrix_dimension_;
    Vector<T> main_diagonal = main_diagonal_;
    Vector<T> gamma = gamma_;

    if (!is_cyclic_) {
        Kokkos::parallel_for(
            "SolveDiagonalNonCyclic",
            Kokkos::RangePolicy<
                Kokkos::DefaultExecutionSpace>(
                0,
                effective_batch_count),
            KOKKOS_LAMBDA(const int k) {
                // ----------------------------------- //
                // Obtain offset for the current batch //
                int batch_idx =
                    batch_stride * k + batch_offset;

                int offset =
                    batch_idx * matrix_dimension;

                // ---------------- //
                // Diagonal Scaling //
                for (int i = 0;
                     i < matrix_dimension;
                     i++) {
                    rhs(offset + i) /=
                        main_diagonal(offset + i);
                }
            });
    }
    else {
        Kokkos::parallel_for(
            "SolveDiagonalCyclic",
            Kokkos::RangePolicy<
                Kokkos::DefaultExecutionSpace>(
                0,
                effective_batch_count),
            KOKKOS_LAMBDA(const int k) {
                // ----------------------------------- //
                // Obtain offset for the current batch //
                int batch_idx =
                    batch_stride * k + batch_offset;

                int offset =
                    batch_idx * matrix_dimension;

                // ---------------- //
                // Diagonal Scaling //
                rhs(offset + 0) /=
                    main_diagonal(offset + 0) +
                    gamma(batch_idx);

                for (int i = 1;
                     i < matrix_dimension;
                     i++) {
                    rhs(offset + i) /=
                        main_diagonal(offset + i);
                }
            });
    }

    Kokkos::fence();
}


private:
int matrix_dimension_;
int batch_count_;
int num_steps_;


Vector<T> main_diagonal_;
Vector<T> sub_diagonal_;

/*
 * Retained for compatibility with the existing solver object.
 *
 * PCR keeps the cyclic auxiliary RHS in team scratch during solve(),
 * so this persistent buffer is no longer required by the algorithm.
 */
Vector<T> buffer_;

Vector<T> gamma_;

/*
 * Strategy A:
 *
 * setup() computes the complete coefficient PCR trajectory once.
 * solve() only replays these multipliers against the RHS.
 */
Vector<T> k1_trajectory_;
Vector<T> k2_trajectory_;

bool is_cyclic_;
bool is_factorized_;

};

} // namespace gmgpolar
